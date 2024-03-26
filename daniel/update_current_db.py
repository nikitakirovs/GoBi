import argparse
import os
import shutil
import sys
import subprocess

class EmptyDirectoryError(Exception):
    '''Exception raised for empty dir. Used in list_subdirectories for -d parameter'''
    def __init__(self, directory):
        self.directory = directory
        super().__init__(f'\nThe directory specified with -d does not contain subdirectories.\n'+
                         f'Please make sure each db to be aliased is contained in its own subdirectory in {directory}.')


def create_fresh_dir(directory_name):
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    else:
        #if input(f'The directory {directory_name} exist already. Its current content will be deleted.'
        #         + '\n If you would like to continue please type yes, y, YES or Y: ') in ['y','yes', 'YES', 'Y']:

            shutil.rmtree(directory_name)
            os.makedirs(directory_name)
        #else:
         #   print('Aborting...')
         #   sys.exit(0)


def list_subdirectories(directory):
    subdirectories = [d for d in os.listdir(directory) if os.path.isdir(os.path.join(directory, d))]
    if not subdirectories:
        raise EmptyDirectoryError(directory)
    return subdirectories

def get_dir_argument(args_dir_dbs):
    count_slashes = args_dir_dbs.count('/')
    if args_dir_dbs.startswith('/'):
        if os.path.exists(args_dir_dbs):
            return args_dir_dbs
        else:
            raise ValueError('Not a valid path.')
    elif '/' in args_dir_dbs:
        if count_slashes == 1 and args_dir_dbs.endswith('/'):
            return args_dir_dbs
        else:
            raise ValueError('Directory name must not contain '/' character unless its the last symbol.')
    else:
        cwd = os.getcwd()
        dir_in_cwd = cwd + f'/{args_dir_dbs}'
        if os.path.exists(dir_in_cwd):
            return dir_in_cwd
        else:
            raise ValueError('Path to directory containing databases does not exist')



def get_dblist(dir_argument):
    if not dir_argument.endswith('/'):
        dir_argument = dir_argument + '/'
    subdirs = list_subdirectories(dir_argument)
    dblist = []
    for dir in subdirs:
        files = os.listdir(dir_argument+f'/{dir}')
        first_file = files[0]
        prefix = first_file.split('.')[0]
        db = dir_argument + f'{dir}/{prefix}'
        dblist.append(db)
    return ' '.join(dblist)

def get_alias_argument(args_aliasdb):
    if args_aliasdb.startswith('/'):
        raise ValueError(f'{args_aliasdb} is not a directory the current working directory')
    elif '/' in args_aliasdb and not args_aliasdb.endswith('/'):
        raise ValueError('Directory name must not contain '/' character unless it\'s the last symbol.')
    else:
        cwd = os.getcwd()
        if args_aliasdb.endswith('/'):
            stripped = args_aliasdb.strip('/')
            return cwd + '/' + stripped
        else:
            return cwd +f'/{args_aliasdb}'

def verify_out_dir(out_dir):
    potential_path = os.path.abspath(out_dir)
    cwd = os.getcwd()

    if potential_path.startswith(cwd):
        return potential_path
    else:
        raise ValueError(f'{potential_path} is not a subdirectory of the current working directory {cwd}')

if __name__ == '__main__':
    ### PARSE
    parser = argparse.ArgumentParser(description=
                                 'Update blast database by creating a new alias db which includes all db in the dir specified with -d')
    parser.add_argument('-d', '--dir_dbs', required=True, help='Path to directory containing one subdirectory for each blastdb')
    parser.add_argument('-a', '--aliasdb', required=True, help='Name of a directory for the alias database in the CURRENT working directory')
    args = parser.parse_args()

    ### SETUP ARGUMENTS FOR blastdb_aliastool

    out_dir = get_alias_argument(args.aliasdb)
    out_dir = verify_out_dir(out_dir)
    out_argument = out_dir + '/aliasdb'
    title_argument = 'alias-db'
    dbtype_argument = 'nucl'
    dir_argument = get_dir_argument(args.dir_dbs)
    dblist_argument = get_dblist(dir_argument)
    
    # blastdb_aliastool -dblist 'db/apis-mellifera bombus/bombus' -dbtype nucl -title current -out current/current
    # both -title and current/<somename> needed

    ### SETUP
    create_fresh_dir(out_dir)
    cmd = (f'blastdb_aliastool -dblist "{dblist_argument}" -dbtype {dbtype_argument} '
          + f'-title {title_argument} -out {out_argument}')
    cmd_output = subprocess.run(cmd, shell=True, capture_output=True)
    print(cmd)
    print(cmd_output.stderr.decode())
