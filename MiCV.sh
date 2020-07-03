#!/bin/bash

## Startup the REDIS server
redis-server &

## Setup the celery task queue
celery multi start worker0 -A tasks.celery:task_queue &

## Start the MiCV web server
gunicorn -b 0.0.0.0:8000 --workers=8 --worker-class gevent --max-requests=10 --timeout=720 index:server
