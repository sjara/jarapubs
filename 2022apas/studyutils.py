"""
Generic methods and classes used throughout.
"""

import os
import pandas as pd
import studyparams
from jaratoolbox import settings


def load_sessions_dataframe(subject):
    """
    Load CSV file containing sessions info.
    """
    sessionsDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, 'sessions')
    sessionsFile = os.path.join(sessionsDir, f'{subject}.csv')
    dframe = pd.read_csv(sessionsFile, index_col=0)
    return dframe
    
def get_sessions(subject, stage):
    """
    Return list of sessions for a given stage.
    """
    dframe = load_sessions_dataframe(subject)
    return list(dframe[dframe.stage==stage].session)


if __name__ == '__main__':
    subject = 'pamo009'; stage = 4
    sessionsList = get_sessions(subject, stage)
    print(sessionsList)
