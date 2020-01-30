cd ~/peptidesim/package
pip install --user -e .
cd ~/scratch
echo "Run this command to use tests:\npython -m pytest -x ../peptidesim/package/tests/"
bash