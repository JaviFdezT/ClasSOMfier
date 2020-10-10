#This script builds the code

#deleted existing directories
rm -rf build
rm -rf CHANGES.txt
rm -rf classifier*
rm -rf classomfier.egg-info/
rm -rf data/
rm -rf dist/
rm -rf nohup.ou

#builds and installs
python setup.py  bdist_wheel;pip install . 

#runs test
cd ../..
#python -c "import classomfier; from classomfier import ClasSOMfier; nn=ClasSOMfier(6.43718,2,'dumpvoid.file',epochs=1);nn.execute()"
#python -c "import classomfier; from classomfier import ClasSOMfier; nn=ClasSOMfier#(6.43718,2,'dumpvoid.file',epochs=1,usenomatrix=True);nn.execute()"
cd -



