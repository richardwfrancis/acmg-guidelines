#!/usr/bin/env python3

import requests
import json

rsids='rs121434224, rs587776569, rs121434225, rs121434226, rs121434227, rs28942093, rs28942094, rs121434228, rs28942095, rs121434229, rs121434230, rs121434231, rs387906967, rs121434227, rs387906968, rs397515434, rs200166664, rs886038202, rs863225063, rs879255573, rs879255243, rs752673677, rs879255244, rs879255245, rs2004640, rs10954213, rs2070197'

for rsid in rsids.split(','):
    rsid = rsid.strip()
    print( rsid, flush=True )
    url = 'https://grch37.rest.ensembl.org/variation/human/'+rsid
    # return type is JSON so bring the data back to an associative array
    r_ensembl = requests.get(url, data={'content-type':'application/json'}, verify=True)
    if ( r_ensembl.status_code == 200 ):
        parsed_json_ensembl = json.loads(r_ensembl.text)
        # interesting fields
        for param in ['mappings','clinical_significance','synonyms']:
            # again check that the field exists in the record
            if param in parsed_json_ensembl:
                # if it's in there then process the contents
                for data in parsed_json_ensembl[param]:
                    # the mappings record contains the location information we want
                    if param == 'mappings':
                        for coord in ['seq_region_name','strand','start','end','allele_string']:
                            print ( "\t{} => {}".format(coord,data[coord]) , flush=True)
                    # the clinical_significance and synonyms records are arrays
                    else:
                        print ( "\t{} => {}".format(param,data), flush=True )
            else:
                print("\t{} not in this record".format(param), flush=True)
    else:
        print("\tfailed with error code {}: {}".format(r_ensembl.status_code, url), flush=True)
        continue
