from __future__ import print_function
import requests
import xml.etree.ElementTree as ET
import numpy

"""
The following lists the arguments that can be used:

mip_era, activity_id, model_cohort, product, source_id, institution_id, source_type, source_id, institution_id, source_type, nominal_resolution, experiment_id, sub_experiment_id, variant_label, grid_label, table_id, frequency, realm, variable_id, cf_standard_name, data_node (for the specific entries refer to the search page directly: https://esgf-node.llnl.gov/search/cmip6/

For example, the following line may be inserted in the code to obtain a list of netcdfs that meet the required options or arguments (as above):

result = esgf_search(activity_id='CMIP', table_id='Amon', variable_id='tas', experiment_id='historical',institution_id="NCAR", source_id="CESM2", member_id="r10i1p1f1")

Script source: http://gallery.pangeo.io/repos/pangeo-gallery/cmip6/search_and_load_with_esgf_opendap.html
Author: Unknown
I got the original version from a word document published by ESGF
https://docs.google.com/document/d/1pxz1Kd3JHfFp8vR2JCVBfApbsHmbUQQstifhGNdc6U0/edit?usp=sharing

API AT: https://github.com/ESGF/esgf.github.io/wiki/ESGF_Search_REST_API#results-pagination
"""


def esgf_search(
    server="https://esgf-node.llnl.gov/esg-search/search",
    files_type="OPENDAP",
    local_node=True,
    project="CMIP6",
    verbose=False,
    format="application%2Fsolr%2Bjson",
    use_csrf=False,
    **search
):
    client = requests.session()
    payload = search
    payload["project"] = project
    payload["type"] = "File"
    if local_node:
        payload["distrib"] = "false"
    if use_csrf:
        client.get(server)
        if "csrftoken" in client.cookies:
            # Django 1.6 and up
            csrftoken = client.cookies["csrftoken"]
        else:
            # older versions
            csrftoken = client.cookies["csrf"]
        payload["csrfmiddlewaretoken"] = csrftoken

    payload["format"] = format

    offset = 0
    numFound = 10000
    all_files = []
    files_type = files_type.upper()
    while offset < numFound:
        payload["offset"] = offset
        url_keys = []
        for k in payload:
            url_keys += ["{}={}".format(k, payload[k])]

        url = "{}/?{}".format(server, "&".join(url_keys))
        print(url)
        r = client.get(url)
        r.raise_for_status()
        resp = r.json()["response"]
        numFound = int(resp["numFound"])
        resp = resp["docs"]
        offset += len(resp)
        for d in resp:
            if verbose:
                for k in d:
                    print("{}: {}".format(k, d[k]))
            url = d["url"]
            for f in d["url"]:
                sp = f.split("|")
                if sp[-1] == files_type:
                    all_files.append(sp[0].split(".html")[0])
    return sorted(all_files)
