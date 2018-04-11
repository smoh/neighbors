from abc import ABCMeta, abstractmethod
import requests

__all__ = [
    ARIADQLValidator,
    CDSADQLValidator
]


class ADQLValidator(object):
    """Base class to share annotate method"""

    @classmethod
    def annotate(cls, query, tables=None):
        """Annotate ADQL query with ADQL validator

        query : str
            query string
        tables : str
            possible options depend on the service

        Prints the annotated query
        """

        errors = cls.check(query, tables=tables)

        if errors is None:
            return

        annotation_str_chars = [[' '] * len(line) for line in query.splitlines()]
        for ie, e in enumerate(errors):
            beginline, begincolumn, endline, endcolumn = \
                e['beginLine']-1, e['beginColumn']-1, e['endLine']-1, e['endColumn']-1
            if beginline == endline:
                annotation_str_chars[beginline][begincolumn:endcolumn+1] = '-'*(endcolumn-begincolumn+1)
            elif beginline < endline:
                for i in range(beginline, endline+1):
                    begin, end = 0, len(annotation_str_chars[i])
                    if i == beginline:
                        begin = begincolumn
                    if i == endline:
                        end = endcolumn
                    annotation_str_chars[i][begin:end+1] = '-'*(end-begin+1)
            annotation_str_chars[beginline][begincolumn] = str(ie)
        annotation_str_lines = [''.join(line) for line in annotation_str_chars]
        annotated_query = []
        for q, a in zip(query.splitlines(), annotation_str_lines):
            annotated_query.append(q)
            if len(a.strip()) != 0:
                annotated_query.append(a)
        errorlines = ['{:2d}: {:s}'.format(ie, e['message']) for ie, e in enumerate(errors)]
        message = '\n'.join(annotated_query + errorlines)
        print(message)


class ARIADQLValidator(ADQLValidator):
    """Validator provided by the Gaia ARI data center"""

    possible_tables = [None, 'ari', 'esa']

    @staticmethod
    def check(query, tables=None):

        if tables not in ARIADQLValidator.possible_tables:
            raise ValueError("tables must be one of {!r}".format(possible_tables))
        query = query.strip()
        payload = dict(query=query, target=tables)
        r = requests.get(
            "http://gaia.ari.uni-heidelberg.de/adql-validator/parse",
            params=payload)

        if r.json()['validation'] == 'ok':
            return
        elif r.json()['validation'] == 'error':
            errors = r.json()['errors']
            return errors


class CDSADQLValidator(ADQLValidator):
    """Validator provided by the CDS"""

    predefined_tables = {
        'gaia': 'http://gea.esac.esa.int/tap-server/tap/tables',
        'gavo': 'http://dc.zah.uni-heidelberg.de/__system__/tap/run/tap/tables',
        'simbad': 'http://simbad.u-strasbg.fr/simbad/sim-tap/tables',
        'vizier': 'http://tapvizier.u-strasbg.fr/TAPVizieR/tap/tables'
    }

    @staticmethod
    def check(query, tables=None):

        tableurl = CDSADQLValidator.predefined_tables[tables] \
            if tables in CDSADQLValidator.predefined_tables \
            else tables
        #NOTE: if query string starts with "\n", I get an error
        #'Connection reset by peer'. Unclear why.
        query = query.strip()
        payload = dict(query=query, tables=tableurl)
        r = requests.get(
            'http://cdsportal.u-strasbg.fr/adqltuto/adqlvalidate',
            params=payload)

        if r.json()['validation'] == 'ok':
            return
        elif r.json()['validation'] == 'error':
            errors = r.json()['errors']
            return errors


def test_validators():

    query = """SELECT COUNT(*) FROM gaiadr1.tgas_source"""
    print(ARIADQLValidator.check(query))
    ARIADQLValidator.annotate(query)
    # assert ARIADQLValidator.check(query) is None

    print(CDSADQLValidator.check(query))
    CDSADQLValidator.annotate(query)

def test_check_adql():
    query = """SELECT COUNT(*) FROM gaiadr1.tgas_source"""
    assert CDSADQLValidator.check(query) is None

    query = """SELECT COUNT(*) FROM gaiadr1.tgas_source"""
    assert CDSADQLValidator.check(query, tables='gaia') is None

    query = "select top 10 *\nfrom gaiadr1.nonexistent\nwhere abc   = d"
    CDSADQLValidator.annotate(query, tables='gaia')
    # expected output:
    #     select top 10 *
    # from gaiadr1.nonexistent
    #      0-------------------
    # where abc   = d
    #       1---    2-
    #  0: Unknown table "gaiadr1.nonexistent" !
    #  1: Unknown column "abc" !
    #  2: Unknown column "d" !

if __name__ == '__main__':
    test_validators()
    test_check_adql()
