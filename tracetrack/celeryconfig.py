import os

broker_url = os.environ.get('CELERY_BROKER_URL', 'redis://redis:6379/0'),
result_backend = os.environ.get('CELERY_RESULT_BACKEND', 'redis://redis:6379/0')

accept_content = ['pickle']
task_serializer = 'pickle'
result_serializer = 'pickle'

redis_backend_health_check_interval = 30

