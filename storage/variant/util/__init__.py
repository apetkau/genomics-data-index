from pathlib import Path
from mimetypes import guess_type
from os.path import basename, splitext


def get_genome_name(file: Path) -> str:
    '''
    Gets the genome name (filename minus extension). Accounts for gzipped/non-gzipped files.
    :param file: The file.
    :return: The genome name from the file.
    '''
    encoding = guess_type(str(file))[1]

    if encoding == 'gzip':
        ref_name = splitext(basename(file).rstrip('.gz'))[0]
    else:
        ref_name = splitext(basename(file))[0]

    return ref_name