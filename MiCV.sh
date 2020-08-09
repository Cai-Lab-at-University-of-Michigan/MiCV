#!/bin/bash

## Startup the REDIS server
redis-server &

## Setup the celery task queue
celery multi start worker0 -A tasks.celery:task_queue &

## Start the MiCV web server
gunicorn -b 0.0.0.0:8000 --workers=4 --threads=2 --worker-class gevent --max-requests=25 --timeout=1440 index:server