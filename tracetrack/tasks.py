from celery import Celery
from tracetrack.alignment_utils import Settings
from typing import List
from tracetrack.entities.alignment import Alignment

celery = Celery('tasks')
config_obj = 'tracetrack.celeryconfig'
celery.config_from_object(config_obj)


class AlignmentTaskResult:
    alignment: Alignment
    warnings: List[str]
    settings: Settings

    def __init__(self, alignment, warnings, settings):
        self.alignment = alignment
        self.warnings = warnings
        self.settings = settings


@celery.task(name='tasks.align_sequences')
def align_sequences_celery(seqlist, ref_id, population_name, db, settings, warnings_orig):
    alignment, warnings = db.align_sequences(seqlist, ref_id, population=population_name)
    warnings += warnings_orig
    return AlignmentTaskResult(alignment=alignment, warnings=warnings, settings=settings)


@celery.task(name='tasks.save_submission')
def save_submission(sequences, population_names, db):
    print("Started")
    for seqlist in sequences:
        for record in seqlist:
            db.match_trace_to_ref(record)
    print("Ready")
    return sequences, population_names, db
