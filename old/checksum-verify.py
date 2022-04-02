# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 19:39:25 2020

@author: Naim
"""

import argparse
import base64


def parse_param():
    parser = argparse.ArgumentParser(description="Fix iRODS sha2 checksum manifest, and/or optionally compare against sha256sum manifest")
    parser.add_argument('ichksum', help='iRODS ichksum output')
    parser.add_argument('-i','--irods_strm', help="String to remove from iRODS ichksum directory output; \
                        skip the initial 'C- ' for this input; must start with '/resarchivezone'", default = '')
    parser.add_argument('-p','--parentdirfilesFilename', help="Filename containing the names of the files in uppermost dir. \
                        One can obtain this by running ls -p |grep -v '/$'")
    parser.add_argument('-m','--manifestFile', help="If supplied, output file will give the differences")
    parser.add_argument('-o','--out', help="output filename for fixed ichksum manifest or \
                        if the manifest file specified via -m option then outputs the missing/non-sha256 matching files")
    args = parser.parse_args()
    ichksum = args.ichksum
    irods_strm = args.irods_strm
    parentdirfilesFilename = args.parentdirfilesFilename
    manifestFile = args.manifestFile
    out = args.out
    
    return (ichksum, irods_strm, parentdirfilesFilename, manifestFile, out)


if __name__=='__main__':
    ichksum, irods_strm, parentdirfilesFilename, manifestFile, out = parse_param()
    
    print(f'ichksum: {ichksum}')
    print(f'irods_strm: {irods_strm}')
    print(f'parentdirfilesFilename: {parentdirfilesFilename}')
    print(f'manifestFile: {manifestFile}')
    print(f'out: {out}')
    
    # ichksum = 'irods-manifest-sha256sum.txt'
    # irods_strm = "/resarchivezone/home/struglis/datasets/epilepsy/BIOJUME/2019-05-EPIC_data_from_Deb_Sophie/"
    # parentdirfilesFilename = 'parentdirfiles.txt'
    # manifestFile = 'manifest-sha256.txt'
    # out = ichksum.replace('.txt', '_fix.txt')
    
    parentdirfiles = []
    if parentdirfilesFilename is not None:
        with open(parentdirfilesFilename, 'r') as f:
            for line in f:
                parentdirfiles.append(line.strip())

    fixed_ichksum = ''
    sha_dict = dict({})
    with open(ichksum, 'r') as f:
        dirstr = './'
        for line in f:
            if line.startswith('C- '):
                dirstr = './' + line.replace('C- ','').replace(irods_strm, '').strip().replace(':','/')
            else:
                linesplit = line.strip().split(' ')
                afile = linesplit[0:(len(linesplit)-2)]
                afile = ' '.join(afile).strip()
                if afile not in parentdirfiles:
                    afile = dirstr + afile
                else:
                    afile = './' + afile
                sha = linesplit[len(linesplit)-1]
                if sha.startswith('sha2:'):
                    sha = sha.replace('sha2:','')
                    sha = base64.b64decode(sha).hex()
                    fixedline = sha + '  ' + afile
                    fixed_ichksum += fixedline + '\r\n'
                    sha_dict[sha] = afile
    
    
    if manifestFile is not None:
        manifest_sha_dict = dict({})
        with open(manifestFile, 'r') as f:
            for line in f:
                linesplit = line.strip().split(' ')
                sha = linesplit[0]
                afile = ' '.join(linesplit[1:len(linesplit)]).strip()
                manifest_sha_dict[sha] = afile
        missing_shas = [ x for x in manifest_sha_dict.keys() if x not in sha_dict.keys() ]
        fixed_ichksum = ""
        if len(missing_shas) > 0:
            for miss_sha in missing_shas:
                fixed_ichksum += miss_sha + '  ' + manifest_sha_dict[miss_sha] + '\r\n'
            

    if out is not None:
        with open(out, 'w') as f:
            f.write(fixed_ichksum)
    else:
        print(fixed_ichksum)