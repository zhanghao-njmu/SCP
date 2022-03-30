"""Class for reading, parsing, and downloading data from the Harmonizome API.
"""

import gzip
import json
import os
import logging

# Support for both Python2.X and 3.X.
# -----------------------------------------------------------------------------
try:
    from io import BytesIO
    from urllib.request import urlopen
    from urllib.error import HTTPError
    from urllib.parse import quote_plus
except ImportError:
    from StringIO import StringIO as BytesIO
    from urllib2 import urlopen, HTTPError
    from urllib import quote_plus

try:
    input_shim = raw_input
except NameError:
    # If `raw_input` throws a `NameError`, the user is using Python 2.X.
    input_shim = input


# Enumerables and constants
# -----------------------------------------------------------------------------

class Enum(set):
    """Simple Enum shim since Python 2.X does not have them.
    """

    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError


# The entity types supported by the Harmonizome API.
class Entity(Enum):

    DATASET = 'dataset'
    GENE = 'gene'
    GENE_SET = 'gene_set'
    ATTRIBUTE = 'attribute'
    GENE_FAMILY = 'gene_family'
    NAMING_AUTHORITY = 'naming_authority'
    PROTEIN = 'protein'
    RESOURCE = 'resource'


def json_from_url(url):
    """Returns API response after decoding and loading JSON.
    """
    response = urlopen(url)
    data = response.read().decode('utf-8')
    return json.loads(data)


VERSION = '1.0'
API_URL = 'https://maayanlab.cloud/Harmonizome/api'
DOWNLOAD_URL = 'https://maayanlab.cloud/static/hdfs/harmonizome/data'

# This config objects pulls the names of the datasets, their directories, and
# the possible downloads from the API. This allows us to add new datasets and
# downloads without breaking this file.
config = json_from_url(API_URL + '/dark/script_config')
DOWNLOADS = [x for x in config.get('downloads')]
DATASET_TO_PATH = config.get('datasets')


# Harmonizome class
# -----------------------------------------------------------------------------

class Harmonizome(object):

    __version__ = VERSION
    DATASETS = DATASET_TO_PATH.keys()

    @classmethod
    def get(cls, entity, name=None, start_at=None):
        """Returns a single entity or a list, depending on if a name is
        provided. If no name is provided and start_at is specified, returns a
        list starting at that cursor position.
        """
        if name:
            name = quote_plus(name)
            return _get_by_name(entity, name)
        if start_at is not None and type(start_at) is int:
            return _get_with_cursor(entity, start_at)
        url = '%s/%s/%s' % (API_URL, VERSION, entity)
        result = json_from_url(url)
        return result

    @classmethod
    def next(cls, response):
        """Returns the next set of entities based on a previous API response.
        """
        start_at = _get_next(response)
        entity = _get_entity(response)
        return cls.get(entity=entity, start_at=start_at)

    @classmethod
    def download(cls, datasets=None, what=None):
        """For each dataset, creates a directory and downloads files into it.
        """
        # Why not check `if not datasets`? Because in principle, a user could 
        # call `download([])`, which should download nothing, not everything.
        # Why might they do this? Imagine that the list of datasets is
        # dynamically generated in another user script.
        if datasets is None:
            datasets = cls.DATASETS
            warning = 'Warning: You are going to download all Harmonizome '\
                      'data. This is roughly 30GB. Do you accept?\n(Y/N) '
            resp = input_shim(warning)
            if resp.lower() != 'y':
                return

        for dataset in datasets:
            if dataset not in cls.DATASETS:
                msg = '"%s" is not a valid dataset name. Check the `DATASETS`'\
                      ' property for a complete list of names.' % dataset
                raise AttributeError(msg)
            if not os.path.exists(dataset):
                os.mkdir(dataset)

            if what is None:
                what = DOWNLOADS

            for dl in what:
                path = DATASET_TO_PATH[dataset]
                url = '%s/%s/%s' % (DOWNLOAD_URL, path, dl)

                try:
                    response = urlopen(url)
                except HTTPError as e:
                    # Not every dataset has all downloads.
                    if what is not None:
                        raise Exception('Error downloading from %s: %s' % (url, e))

                filename = '%s/%s' % (dataset, dl)
                filename = filename.replace('.gz', '')

                if response.code != 200:
                    raise Exception('This should not happen')
                
                if os.path.isfile(filename):
                    logging.info('Using cached `%s`' % (filename))
                else:
                    logging.info('Downloading `%s`' % (filename))
                    _download_and_decompress_file(response, filename)

                yield filename

    @classmethod
    def download_df(cls, datasets=None, what=None, sparse=False, **kwargs):
        for file in cls.download(datasets, what):
            if sparse:
                yield _read_as_sparse_dataframe(file, **kwargs)
            else:
                yield _read_as_dataframe(file, **kwargs)

# Utility functions
# -------------------------------------------------------------------------

def _get_with_cursor(entity, start_at):
    """Returns a list of entities based on cursor position.
    """
    url = '%s/%s/%s?cursor=%s' % (API_URL, VERSION, entity,str(start_at))
    result = json_from_url(url)
    return result


def _get_by_name(entity, name):
    """Returns a single entity based on name.
    """
    url = '%s/%s/%s/%s' % (API_URL, VERSION, entity, name)
    return json_from_url(url)


def _get_entity(response):
    """Returns the entity from an API response.
    """
    path = response['next'].split('?')[0]
    return path.split('/')[3]


def _get_next(response):
    """Returns the next property from an API response.
    """
    if response['next']:
        return int(response['next'].split('=')[1])
    return None


# This function was adopted from here: http://stackoverflow.com/a/15353312.
def _download_and_decompress_file(response, filename):
    """Downloads and decompresses a single file from a response object.
    """
    compressed_file = BytesIO(response.read())
    decompressed_file = gzip.GzipFile(fileobj=compressed_file)

    with open(filename, 'wb+') as outfile:
        outfile.write(decompressed_file.read())


def _getfshape(fn, row_sep='\n', col_sep='\t', open_args={}):
    ''' Fast and efficient way of finding row/col height of file '''
    with open(fn, 'r', newline=row_sep, **open_args) as f:
        col_size = f.readline().count(col_sep) + 1
        row_size = sum(1 for line in f) + 1
        return (row_size, col_size)

def _parse(fn, column_size=3, index_size=3, shape=None,
          index_fmt=None, data_fmt=None,
          index_dtype=None, data_dtype=None,
          col_sep='\t', row_sep='\n',
          open_args={}):
    '''
    Smart(er) parser for processing matrix formats. Evaluate size and construct
     ndframes with the right size before parsing, this allows for more efficient
     loading of sparse dataframes as well. To obtain a sparse representation use:
         data_fmt=scipy.lil_matrix
    This only works if all of the data is of the same type, if it isn't a float
     use:
         data_dtype=np.float64
    
    Returns:
        (column_names, columns, index_names, index, data)
    '''
    import numpy as np

    if index_fmt is None: index_fmt = np.ndarray
    if data_fmt is None: data_fmt = np.ndarray
    if index_dtype is None: index_dtype = np.object
    if data_dtype is None: data_dtype = np.float64

    if shape is not None:
        rows, cols = shape
    else:
        rows, cols = _getfshape(fn, row_sep=row_sep, col_sep=col_sep, open_args=open_args)

    columns = index_fmt((column_size, cols - index_size), dtype=index_dtype)
    index = index_fmt((rows - column_size, index_size), dtype=index_dtype)
    data = data_fmt((rows - column_size, cols - index_size), dtype=data_dtype)

    with open(fn, 'r', newline=row_sep, **open_args) as fh:
        header = np.array([next(fh).strip().split(col_sep)
                           for _ in range(column_size)])

        column_names = header[:column_size, index_size - 1]
        index_names = header[column_size - 1, :index_size]

        columns[:, :] = header[:column_size, index_size:]

        for ind, line in enumerate(fh):
            lh = line.strip().split(col_sep)
            index[ind, :] = lh[:index_size]
            data[ind, :] = lh[index_size:]

        return (column_names, columns, index_names, index, data)

def _parse_df(fn, sparse=False, default_fill_value=None,
             column_apply=None, index_apply=None, df_args={},
             **kwargs):
    import numpy as np
    import pandas as pd
    from scipy.sparse import lil_matrix

    data_fmt = lil_matrix if sparse else np.ndarray
    df_type = pd.SparseDataFrame if sparse else pd.DataFrame
    (
        column_names, columns,
        index_names, index,
        data,
    ) = _parse(fn, data_fmt=data_fmt, **kwargs)

    if column_apply is not None:
        column_names, columns = column_apply(column_names.T, columns.T)
    else:
        column_names, columns = (column_names.T, columns.T)

    if index_apply is not None:
        index_names, index = index_apply(index_names, index)

    return df_type(
        data=data.tocsr() if sparse else data,
        index=pd.Index(
            data=index,
            name=str(index_names),
            dtype=np.object,
        ),
        columns=pd.Index(
            data=columns,
            name=str(column_names),
            dtype=np.object,
        ),
        **df_args,
    )

def _df_column_uniquify(df):
    df_columns = df.columns
    new_columns = []
    for item in df_columns:
            counter = 0
            newitem = item
            while newitem in new_columns:
                    counter += 1
                    newitem = "{}_{}".format(item, counter)
            new_columns.append(newitem)
    df.columns = new_columns
    return df

def _json_ind_no_slash(ind_names, ind):
    return (
        json.dumps([ind_name.replace('/', '|')
                    for ind_name in ind_names]),
        [json.dumps([ii.replace('/', '|')
                     for ii in i])
         for i in ind],
    )

def _read_as_dataframe(fn):
    ''' Standard loading of dataframe '''
    if fn.endswith('gene_attribute_matrix.txt'):
        return _df_column_uniquify(_parse_df(
            fn,
            sparse=False,
            index_apply=_json_ind_no_slash,
            column_apply=_json_ind_no_slash,
            open_args=dict(encoding="latin-1"),
        ))
    elif fn.endswith('gene_list_terms.txt') or fn.endswith('attribute_list_entries.txt'):
        import pandas as pd
        return pd.read_table(fn, encoding="latin-1", index_col=None)
    else:
        raise Exception('Unable to parse this file into a dataframe.')

def _read_as_sparse_dataframe(fn, blocksize=10e6, fill_value=0):
    ''' Efficient loading sparse dataframe '''
    if fn.endswith('gene_attribute_matrix.txt'):
        return _df_column_uniquify(_parse_df(
            fn,
            sparse=True,
            index_apply=_json_ind_no_slash,
            column_apply=_json_ind_no_slash,
            df_args=dict(default_fill_value=0),
            open_args=dict(encoding="latin-1"),
        ))
    else:
        raise Exception('Unable to parse this file into a dataframe.')
