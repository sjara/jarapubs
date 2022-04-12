"""
This file contains the default names and queries used for the study on
Passive exposure effects on learning (active vs active+passive)
"""

STUDY_NAME = '2022apas'


# -- List of animals used in this study --

MICE_ACTIVE_ONLY = ['pamo009', 'pamo011', 'pamo012', 'pamo013', 'pamo016',
                    'pamo019', 'pamo022', 'pamo026', 'pamo020']

MICE_ACTIVE_PASSIVE = ['pamo010', 'pamo014', 'pamo015', 'pamo017', 'pamo018',
                       'pamo021', 'pamo023', 'pamo025', 'pamo024']

MICE_ALL = MICE_ACTIVE_ONLY + MICE_ACTIVE_PASSIVE
MICE_ALL.sort()

SESSIONS_RANGE = ['2022-01-14', '2022-03-07']

# -- See reason in the wiki page for this study --
EXCLUDE_SESSIONS = {'pamo022':['2022-03-13'],
                    'pamo023':['2022-01-24'],
                    'pamo026':['2022-01-21']}


