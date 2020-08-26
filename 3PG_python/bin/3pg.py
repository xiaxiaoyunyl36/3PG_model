'''
File: 3pg.py
Author: joeyzhou1984@gmail.com
Description: the entry point of 3-PG model
Created: 2012-01-02
LastModified: 2012-01-02
'''

import sys
sys.path.append('../lib')
from Model3PG import Model3PG


def read_control_file(fpath):
    model = Model3PG(fpath)
    model.initialize()
    return model


def run_3pg(fpath_control):
    # read the control file
    try:
        model_3pg = read_control_file(fpath_control)
    except:
        print('%s is not a valid control file.' % fpath_control)
        return
    model_3pg.run()


def main(argv):
    if len(argv) == 1:
        print('Please provide at least one control file for running the model')
        # run_3pg('run1.yaml')
    else:
        for fpath in argv[1:]:
            run_3pg(fpath)


if __name__ == '__main__':
    main(sys.argv)
