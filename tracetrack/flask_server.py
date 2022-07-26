from flask import Flask, render_template, request, flash, redirect, send_file
from tracetrack.entities.alignment import Alignment
from tracetrack.entities.reference_db import ReferenceDb
from tracetrack.server_utils import unzip_and_get_sequences, number_prefix, schedule_tasks
import os
import time
from tracetrack.server_utils import FileError, get_trace_files, create_spreadsheet, create_table, \
    ResultsError, get_ref_input, is_best_for_reference, get_color_bin, trace_record_dict, position_dict, hash_trace_name
import hashlib
from tracetrack.tasks import celery, save_submission, AlignmentTaskResult
from tracetrack.entities.scheduler import scheduler, use_scheduler
import numpy as np
from typing import List

THRESHOLD = 0
END_THRESHOLD = 0

use_scheduler('celery')


def md5(s):
    return hashlib.md5(s.encode()).hexdigest()


app = Flask(__name__)
app.secret_key = os.urandom(27)
app.jinja_env.globals.update(zip=zip)
app.jinja_env.globals.update(round=round)
app.jinja_env.globals.update(enumerate=enumerate)
app.jinja_env.globals.update(md5=md5)
app.jinja_env.globals.update(is_best_for_reference=is_best_for_reference)
app.jinja_env.globals.update(get_color_bin=get_color_bin)


@app.template_filter('autoversion')
def autoversion_filter(filename):
    # determining fullpath might be project specific
    bin_dir = os.path.dirname(os.path.realpath(__file__))
    fullpath = os.path.join(bin_dir, filename.lstrip('/'))
    try:
        timestamp = str(os.path.getmtime(fullpath))
    except OSError:
        return filename
    newfilename = "{0}?v={1}".format(filename, timestamp)
    return newfilename


@app.route("/")
def intro_get():
    return render_template('welcome.html')


@app.route("/input")
def input_get():
    return render_template(
        'index.html',
        results=None
    )


@app.route("/help")
def help_get():
    return render_template('help.html')


@app.route("/excel/<task_id>", methods=['GET'])
def output_spreadsheet(task_id):
    results: List[AlignmentTaskResult] = scheduler.get_results(task_id, timeout=10)

    errors = [str(r) for r in results if isinstance(r, Exception)]
    results = [r for r in results if isinstance(r, AlignmentTaskResult)]

    alignments = [res.alignment for res in results]
    warnings = errors + [warn for res in results for warn in res.warnings]

    for warning in warnings:
        flash(warning)

    try:
        spreadsheet = create_spreadsheet(alignments)
    except ResultsError as error:
        return render_template('results.html', display=error.message)

    # finally return the file
    timestr = time.strftime("%Y-%m-%d_%H-%M-%S")
    return send_file(spreadsheet, attachment_filename="TraceTrack{}.xlsx".format(timestr), as_attachment=True,
                     cache_timeout=0)


@app.route("/excel/<task_id>/<al_num>", methods=['GET'])
def export_alignment(task_id, al_num):
    results: List[AlignmentTaskResult] = scheduler.get_results(task_id, timeout=10)

    errors = [str(r) for r in results if isinstance(r, Exception)]
    results = [r for r in results if isinstance(r, AlignmentTaskResult)]

    alignments = [res.alignment for res in results]
    warnings = errors + [warn for res in results for warn in res.warnings]

    for warning in warnings:
        flash(warning)
    try:
        table = create_table(alignments)
        spreadsheet = create_spreadsheet([alignments[int(al_num)]], table, int(al_num))
    except ResultsError as error:
        return render_template('results.html', display=error.message)
    # finally return the file
    timestr = time.strftime("%Y-%m-%d_%H-%M-%S")
    return send_file(spreadsheet, attachment_filename="TraceTrack{}.xlsx".format(timestr), as_attachment=True,
                     cache_timeout=0)


@app.route("/trace/<task_id>/<int:alignment_index>", methods=['GET'])
def trace_get(task_id, alignment_index):
    results: List[AlignmentTaskResult] = scheduler.get_results(task_id, timeout=10)
    alignments = [res.alignment for res in results]

    alignment: Alignment = alignments[alignment_index]

    output_dict = {
        'sequences': [trace_record_dict(t) for t in alignment.aligned_traces],
        'reference': str(alignment.aligned_reference.record.seq),
        'mutated': "".join([pos.result for pos in alignment.positions]),
        'translation': alignment.get_translation(),
        'positions': [position_dict(alignment, p) for p in range(len(alignment.positions))],
        'refid': alignment.ref_id
    }
    return output_dict


@app.route("/input", methods=['POST'])
def input_post():
    # loads input files to a temporary directory saved in variable tmp
    example = 'exampleButton' in request.form

    # check if the post request has the file part
    if (not example and 'tracefile1[]' not in request.files) or 'reference' not in request.files:
        flash('No file part')
        print('No file part')
        return redirect(request.url)

    if example:
        reffile = os.path.join(app.root_path, "../data/example/example_sheet.xlsx")
        db = ReferenceDb.read_file(reffile)

        example_tracefiles = os.path.join(app.root_path, "../data/example/traces")
        seq_lists = [unzip_and_get_sequences(example_tracefiles)]
        population_names = ['Example']
    else:
        try:
            refdir, reffile = get_ref_input(request.files["reference"])
        except FileError as error:
            flash(error.message)
            return redirect(request.url)

        db = ReferenceDb.read_file(os.path.join(refdir.name, reffile))
        if len(db.sequences) == 0:
            flash(
                "The reference sheet contains no sequences. Please upload a reference sheet with a header and one or more lines.")
            return redirect(request.url)

        index = 1
        population_names = []
        seq_lists = []
        while True:
            files = request.files.getlist('tracefile{}[]'.format(index))
            if not files or not any(files):
                break
            try:
                tmpdir = get_trace_files(files, index)
            except FileError as error:
                flash(error.message)
                return redirect(request.url)
            name = request.form['population{}'.format(index)]
            index += 1
            # get list of lists of TraceSeqRecord objects
            seq_lists.append(unzip_and_get_sequences(tmpdir.name))
            tmpdir.cleanup()
            population_names.append(name)

    task_id = scheduler.schedule_task(save_submission, seq_lists, population_names, db)

    return redirect("/settings/" + task_id)


@app.route("/settings/<task_id>", methods=['GET'])
@app.route("/settings/<task_id>/<settings_task>", methods=['GET'])
def alignment_settings_get(task_id, settings_task=None):
    if not scheduler.is_result_ready(task_id):
        return render_template("loading.html", text="Your files are being processed. This page refreshes automatically and will display the settings page once it's ready.")

    sequences, populations, db = scheduler.get_result(task_id)
    if sum([len(seq) for seq in sequences]) == 0:
        flash("No valid trace files were uploaded.")
        return redirect("/")
    quals = []
    for seqlist in sequences:
        quals = quals + [np.array(rec.quality) for rec in seqlist]
    quals = np.concatenate(quals)
    bin_size = 2

    histogram, edges = np.histogram(quals, 100 // bin_size, range=(0, 100))
    hist = []
    for i in range(100 // bin_size):
        hist.append({'bin': edges[i], 'value': float(histogram[i])})

    ref_assignment = {}
    other_refs = {}
    for i, seqlist in enumerate(sequences):
        pop = populations[i]
        for j, seq in enumerate(seqlist):
            dir = "Rev" if seq.reverse else "Fwd"
            seq_id = hash_trace_name(i, j)
            if seq.id in ref_assignment:
                ref_assignment[seq_id] = {
                    "name": f"{seq.id} ({pop})", "ref": seq.reference, "pop": pop,
                    "r_f": dir
                }
                other_refs[seq_id] = [ref.id for ref in db.sequences if ref.id != seq.reference]
            else:
                ref_assignment[seq_id] = {"name": seq.id, "ref": seq.reference, "pop": pop, "r_f": dir}
                other_refs[seq_id] = [ref.id for ref in db.sequences if ref.id != seq.reference]
    ref_assignment = dict(sorted(ref_assignment.items()))

    if settings_task:
        results: List[AlignmentTaskResult] = scheduler.get_results(task_id, timeout=10)
        settings = results[0].settings if len(results) else None
        return render_template("settings.html", threshold=settings.threshold, end_threshold=settings.end_threshold,
                               separate=settings.separate, task_id=task_id,
                               histogram=hist, ref_assignment=ref_assignment, references=other_refs)
    return render_template("settings.html", threshold=THRESHOLD, end_threshold=END_THRESHOLD,
                           task_id=task_id, histogram=hist, ref_assignment=ref_assignment, references=other_refs)


@app.route("/settings/<task_id>", methods=['POST'])
@app.route("/settings/<task_id>/<settings_task>", methods=['POST'])
def alignment_settings_post(task_id, settings_task=None):
    threshold = int(request.form['threshold'])
    end_threshold = int(request.form['end_threshold'])
    separate = 'separate' in request.form

    result = scheduler.get_result(task_id)
    if isinstance(result, Exception):
        error = str(result)
        traceback = result.traceback
        return render_template('results.html', error=error, alignments=None, end_threshold=END_THRESHOLD,
                               is_multiple_populations=False, traceback=traceback,
                               threshold=THRESHOLD, task_id=task_id)

    sequences, population_names, db = result
    for i, seqlist in enumerate(sequences):
        for j, seq in enumerate(seqlist):
            seq_id = hash_trace_name(i, j)
            ref = str(request.form[f'ref_for_{seq_id}'])
            seq.assign_reference(ref)
            direction = str(request.form[f'dir_for_{seq_id}'])
            dir_flag = True if direction == "Rev" else False
            seq.flag_as_reverse(dir_flag)

    new_id = schedule_tasks(sequences, population_names, db, separate, threshold, end_threshold)
    return redirect("/results/" + task_id + "/" + new_id)


@app.route("/results/<files_id>/<task_id>", methods=['GET'])
def results_page(files_id, task_id):
    if not scheduler.are_results_ready(task_id):
        num_running, num_completed, total = scheduler.get_results_progress(task_id)
        return render_template("loading.html",
                               text=f"{num_completed} out of {total} alignments are completed. This page refreshes automatically and will display the result once it's ready.")

    results: List[AlignmentTaskResult] = scheduler.get_results(task_id, timeout=10)
    errors = [str(r) for r in results if isinstance(r, Exception)]
    results = [r for r in results if isinstance(r, AlignmentTaskResult)]

    alignments = [res.alignment for res in results]
    warnings = errors + [warn for res in results for warn in res.warnings]
    settings = results[0].settings if len(results) else None

    for warning in warnings:
        flash(warning)
    populations = [al.population for al in alignments]
    pop_number = len(set(populations))

    table = create_table(alignments)

    return render_template('results.html', error=None, alignments=alignments, settings=settings,
                           is_multiple_populations=pop_number, results=table, task_id=task_id, files_id=files_id)
