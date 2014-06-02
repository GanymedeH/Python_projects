from distutils.core import setup
import py2exe


if __name__=="__main__":
    setup(console=['testADP.py','dataReports.py','exeADP.py','networkGenerator.py','simulationEngine.py'])