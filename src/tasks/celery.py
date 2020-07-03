from celery import Celery

task_queue = Celery('MiCV',
                    broker='redis://127.0.0.1:6379',
                    backend='redis://127.0.0.1:6379',
                    include=['tasks.tasks'])
# Optional configuration, see the application user guide.
task_queue.conf.update(
    result_expires=3600,
)
if __name__ == '__main__':
    task_queue.start()