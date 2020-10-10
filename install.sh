rm -rf build
rm -rf CHANGES.txt
rm -rf classifier*
rm -rf classomfier.egg-info/
rm -rf data/
rm -rf dist/
rm -rf nohup.ou

rm -r /home/javi/system/anaconda/lib/python3.8/site-packages/classomf*
python setup.py build bdist_wheel;pip install . 


ls -l /home/javi/system/anaconda/lib/python3.8/site-packages/classomfier

cd classomfier/test
python -c "import classomfier; from classomfier import ClasSOMfier; nn=ClasSOMfier(6.43718,2,'dumpvoid.file');nn.execute()"
cd -



