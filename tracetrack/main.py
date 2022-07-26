import click
import traceback
from tracetrack.web import web


class MainGroup(click.Group):
    def __call__(self, *args, **kwargs):
        try:
            return self.main(*args, **kwargs)
        except Exception as e:
            message = e.args[0] if e.args else ''
            traceback.print_exc()
            click.echo(f'SuprSeqr failed with {type(e).__name__}' + (f': {message}' if message else ''))
            if len(e.args) > 1:
                for arg in e.args[1:]:
                    click.echo(arg, err=True)
            exit(1)


@click.group(cls=MainGroup)
def main():
    pass


# Register all commands
main.add_command(web)


if __name__ == '__main__':
    main()
