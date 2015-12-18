#!/usr/bin/env python

import sys
import os
import re
import tempfile
import subprocess
import fileinput
import shutil
import optparse
import urllib2
from ftplib import FTP
import tarfile

from galaxy.util.json import from_json_string, to_json_string


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)

def get_config_dict(config,dataset_directory=None):
    keys = ['dataset_directory','ensembl_organism','ensembl_prefix','ensembl_version','ensembl_genome_version','ucsc_genome_version','ncbi_organism','ncbi_prefix','chromosomes','mt_chromosome','gene_sources','ig_gene_sources','rrna_gene_sources']
    pat = '^([^=]+?)\s*=\s*(.*)$'
    config_dict = {}
    try:
        fh = open(config)
        for i,l in enumerate(fh):
           line = l.strip() 
           if line.startswith('#'):
               continue
           m = re.match(pat,line)
           if m and len(m.groups()) == 2:
               (k,v) = m.groups()
               if k in keys:
                   config_dict[k] = v
    except Exception, e:
        stop_err( 'Error parsing %s %s\n' % (config,str( e )) )
    else:
        fh.close()
    if dataset_directory:
        config_dict['dataset_directory'] = dataset_directory
    return config_dict

def run_defuse_script(data_manager_dict, params, target_directory, dbkey, description, config, script):
    if not os.path.isdir(target_directory):
        os.makedirs(target_directory)
    ## Name the config consistently with data_manager_conf.xml
    #  copy the config file to the target_directory
    #  when DataManager moves files to there tool-data location, the config will get moved as well,
    #   and the value_translation in data_manager_conf.xml will tell us the new location
    #  defuse.xml will use the path to this config file to set the dataset_directory
    config_name = '%s.config' % dbkey
    defuse_config = os.path.join( target_directory, config_name)
    shutil.copyfile(config,defuse_config) 
    cmd = "/bin/bash %s %s" % (script,target_directory)
    # Run
    try:
        tmp_out = tempfile.NamedTemporaryFile().name
        tmp_stdout = open( tmp_out, 'wb' )
        tmp_err = tempfile.NamedTemporaryFile().name
        tmp_stderr = open( tmp_err, 'wb' )
        proc = subprocess.Popen( args=cmd, shell=True, cwd=".", stdout=tmp_stdout, stderr=tmp_stderr )
        returncode = proc.wait()
        tmp_stderr.close()
        # get stderr, allowing for case where it's very large
        tmp_stderr = open( tmp_err, 'rb' )
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += tmp_stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stdout.close()
        tmp_stderr.close()
        if returncode != 0:
            raise Exception, stderr

        # TODO: look for errors in program output.
    except Exception, e:
        stop_err( 'Error creating defuse reference:\n' + str( e ) )
    config_dict = get_config_dict(config, dataset_directory=target_directory)
    data_table_entry = dict(value=dbkey, dbkey=dbkey, name=description, path=config_name)
    _add_data_table_entry( data_manager_dict, data_table_entry )
def _add_data_table_entry( data_manager_dict, data_table_entry ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables']['defuse_reference'] = data_manager_dict['data_tables'].get( 'defuse_reference', [] )
    data_manager_dict['data_tables']['defuse_reference'].append( data_table_entry )
    return data_manager_dict

def main():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-k', '--dbkey', dest='dbkey', action='store', type="string", default=None, help='dbkey' )
    parser.add_option( '-d', '--description', dest='description', action='store', type="string", default=None, help='description' )
    parser.add_option( '-c', '--defuse_config', dest='defuse_config', action='store', type="string", default=None, help='defuse_config' )
    parser.add_option( '-s', '--defuse_script', dest='defuse_script', action='store', type="string", default=None, help='defuse_script' )
    (options, args) = parser.parse_args()

    filename = args[0]

    params = from_json_string( open( filename ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory )
    data_manager_dict = {}

     
    #Create Defuse Reference Data
    run_defuse_script( data_manager_dict, params, target_directory, options.dbkey, options.description,options.defuse_config,options.defuse_script)

    #save info to json file
    open( filename, 'wb' ).write( to_json_string( data_manager_dict ) )

if __name__ == "__main__": main()

