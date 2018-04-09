import requests

predefined_tables = {
    'gaia': 'http://gea.esac.esa.int/tap-server/tap/tables',
    'gavo': 'http://dc.zah.uni-heidelberg.de/__system__/tap/run/tap/tables',
    'simbad': 'http://simbad.u-strasbg.fr/simbad/sim-tap/tables',
    'vizier': 'http://tapvizier.u-strasbg.fr/TAPVizieR/tap/tables'
}


def check_adql(query, tables=None):
    """Check ADQL query with ADQL validator

    query : str
        query string
    tables : str
        url to TAP servies tables list.
        Can be one of 'gaia', 'gavo', 'simbad', 'vizier'

    Returns `None` if valid else a list or errors (dictionaries).
    """
    tableurl = predefined_tables[tables] if tables in predefined_tables\
        else tables
    payload = dict(query=query, tables=tableurl)
    r = requests.get(
        'http://cdsportal.u-strasbg.fr/adqltuto/adqlvalidate',
        params=payload)
    if r.json()['validation'] == 'ok':
        return
    elif r.json()['validation'] == 'error':
        return r.json()['errors']


def test_check_adql():
    query = """SELECT COUNT(*) FROM gaiadr1.tgas_source"""
    assert check_adql(query) is None

    query = """SELECT COUNT(*) FROM gaiadr1.tgas_source"""
    assert check_adql(query, tables='gaia') is None
