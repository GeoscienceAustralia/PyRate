echo "startup script started"

export PYRATEPATH=$HOME"/PyRate"
export PYTHONPATH=$PYRATEPATH/:$PYTHONPATH

export http_proxy="http://<username>:<password>@proxy.ga.gov.au:8080"
export https_proxy="http://<username>:<password>@proxy.ga.gov.au:8080"