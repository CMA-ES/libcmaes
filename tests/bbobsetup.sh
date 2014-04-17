#!/bin/bash
/usr/bin/wget http://coco.lri.fr/downloads/download13.09/bbobc.tar.gz
tar xvzf bbobc.tar.gz
/usr/bin/python bbob.v13.09/c/createfolders.py cmaes_bbob
/usr/bin/python bbob.v13.09/c/createfolders.py ipop_bbob
/usr/bin/python bbob.v13.09/c/createfolders.py bipop_bbob
/usr/bin/python bbob.v13.09/c/createfolders.py acmaes_bbob
/usr/bin/python bbob.v13.09/c/createfolders.py aipop_bbob
/usr/bin/python bbob.v13.09/c/createfolders.py abipop_bbob
cd ../
./configure --enable-bbob
make
