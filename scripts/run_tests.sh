#!/bin/bash
source ~/PyRateVenv/bin/activate
pip install -U pytest
cd ~/PyRate
pytest tests/
