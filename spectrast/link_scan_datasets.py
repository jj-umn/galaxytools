#!/usr/bin/env python

from __future__ import print_function

import argparse
import os
import os.path
import xml.sax
from difflib import SequenceMatcher


def __main__():
    parser = argparse.ArgumentParser(
        description='link spectrum datasets to the name used' +
                    ' in the identification dataset')
    parser.add_argument(
        'ident_files', nargs='+',
        help='Pepxml or mzIdentML')
    parser.add_argument(
        '-n', '--scan_name', default=[], action='append',
        help='Name for scan file')
    parser.add_argument(
        '-f', '--scan_file', default=[], action='append',
        help='Path for scan file')
    args = parser.parse_args()

    class MzidHandler(xml.sax.ContentHandler):

        def __init__(self):
            xml.sax.ContentHandler.__init__(self)
            self.spectraDataFiles = []
            self.spectraDataNames = []
            self.searchDatabaseFiles = []
            self.searchDatabaseNames = []

        def startElement(self, tag, attrs):
            if tag == 'SpectraData':
                id = attrs['id']
                path = attrs['location']
                filename = os.path.basename(path)
                name = attrs['name'] if 'name' in attrs else None
                self.spectraDataFiles.append(filename)
                self.spectraDataNames.append(name if name else id)
                print ("SpectraData: %s  %s" % (name if name else id, path))
            if tag == 'SearchDatabase':
                id = attrs['id']
                path = attrs['location']
                filename = os.path.basename(path)
                name = attrs['name'] if 'name' in attrs else None
                self.searchDatabaseFiles.append(filename)
                self.searchDatabaseNames.append(name if name else id)
                print ("SearchDatabase: %s  %s" % (name if name else id, path))

        def endElement(self, name):
            pass

        def characters(self, data):
            pass

    class PepXmlHandler(xml.sax.ContentHandler):

        def __init__(self):
            xml.sax.ContentHandler.__init__(self)
            self.spectraDataFiles = []
            self.spectraDataNames = []

        def startElement(self, tag, attrs):
            if tag == 'msms_run_summary':
                basename = attrs['base_name']
                name = os.path.basename(basename)
                ext = attrs['raw_data']
                path = '%s%s' % (basename, ext)
                filename = os.path.basename(path)
                self.spectraDataFiles.append(filename)
                self.spectraDataNames.append(name)
                print ("SpectraData: %s  %s" % (name, path))

        def endElement(self, name):
            pass

        def characters(self, data):
            pass

    parser = xml.sax.make_parser()
    parser.setFeature(xml.sax.handler.feature_namespaces, 0)
    handler = PepXmlHandler()
    parser.setContentHandler(handler)
    for ident in args.ident_files:
        parser.parse(ident)

    spectra_names = handler.spectraDataFiles

    def best_match(name):
        if name in spectra_names:
            return name
        try:
            r = [SequenceMatcher(None, name, spectra_names[x]).ratio()
                 for x in range(len(spectra_names))]
            return spectra_names[r.index(max(r))]
        except Exception, e:
            print ("best_match: %s  %s" % (name, e))

    for i, name in enumerate(args.scan_name):
        path = args.scan_file[i] if len(args.scan_file) > i else ''
        (root, ext) = os.path.splitext(name)
        print ("SpectraFile: %s  %s" % (name, path))
        iname = best_match(name)
        print ("IdentName: %s  %s" % (name, iname))
        if not os.path.exists(iname) and os.path.exists(path):
            os.symlink(path, iname)
            print ("%s -> %s" % (iname, path))


if __name__ == "__main__":
    __main__()
