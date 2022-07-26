import click
from tracetrack.entities.scheduler import use_scheduler
from tracetrack.flask_server import app
import time


@click.command()
@click.option('--host', default='localhost', help='Server port')
@click.option('--port', type=int, default=5001, help='Server port')
def web(host, port):
    """Run the TraceTrack web interface in simplified mode (without a processing queue)"""
    click.echo('Note! This is a simplified TraceTrack server!')
    click.echo('      See README on how to run TraceTrack with a processing queue.')
    click.echo('')
    time.sleep(1)
    click.echo('Starting TraceTrack...')

    use_scheduler('simple')

    app.run(debug=True, host=host, port=port, processes=1, use_reloader=False)
