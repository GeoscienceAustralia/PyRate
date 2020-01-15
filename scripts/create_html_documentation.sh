#!/bin/bash
cd ~/PyRate/docs
make html
echo "Use following cmd to copy file to local machine:"
echo "scp -r `whoami`@gadi.nci.org.au:`pwd ~`/PyRate/docs/_build/html C:/preview"
