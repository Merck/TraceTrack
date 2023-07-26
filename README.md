# TraceTrack

Sanger sequencing analysis web app. View the [publication](https://biorxiv.org/cgi/content/short/2022.07.28.501824v1) or try the app [online](https://tracetrack.dichlab.org).

## Overview
The TraceTrack application serves for automatic batch analysis of Sanger sequencing trace files. It
accepts a set of trace files and a reference file, performs alignments of trace sequences to 
corresponding reference sequences and displays an overview of resulting alignments.

For more information on how to use the app, see the "Help" page in the application itself.


## Run TraceTrack in Docker container

### 1. Install Docker

See https://docs.docker.com/docker-for-mac/install/

### 2. Build all images using Docker Compose

```bash
make docker-build
```

### 3. Run all services using Docker Compose

```bash
make docker-run
```

To build and run, you can use:
```bash
make docker-build docker-run
```
### 4. Handle code updates

After your code is updated, you will need to stop the services, run build and start again. 
See the next section for info on running locally with flask auto-reload.

### 5. Frontend

Frontend is available through your browser at http://localhost:5001/

## Run TraceTrack using Conda

### 1. Install Conda

See https://docs.conda.io/en/latest/miniconda.html

### 2. Install Redis server

You can [install Redis using Brew](https://medium.com/@petehouston/install-and-config-redis-on-mac-os-x-via-homebrew-eb8df9a4f298).

### 3. Setup environment

```bash
# Install dependencies using environment.yml
make install
# Update environment with:
make update-env
```

### 4. Run all services

You will have to run each service in a separate terminal (Use Cmd+T to open a new tab):

```bash
# Run Redis server
make redis

# In a separate terminal, run celery worker queue
make celery

# In a separate terminal, run flask web server
make web
```
### 5. Handle code updates

After your code is updated, the flask web service should refresh automatically. However, the celery service
needs to be stopped and started manually, so you will need to do that if you update code 
that is executed from the workers.

### 6. Frontend

Frontend is available through your browser at http://localhost:5001/

## Files

- [Dockerfile](Dockerfile) 
- [Makefile](Makefile) 
- [README.md](README.md) 
- [tracetrack](tracetrack) 
  - [flask_server.py](tracetrack/flask_server.py) Basic app functionality, processing uploaded files
  - [server_utils.py](tracetrack/server_utils.py) Functionality for displaying, downloading and exporting data
  - [alignment_utils.py](tracetrack/alignment_utils.py) Contains functions for sequence processing and alignment using Clustal Omega
  - [tasks.py](tracetrack/tasks.py) Task for task queue
  - [celeryconfig.py](tracetrack/celeryconfig.py) Configuration of task queue
  - [loading_data.py](tracetrack/loading_data.py) Loading and processing input files functionality
  - [web.py](tracetrack/web.py) CLI for running with no task queue
  - [entities](tracetrack/entities)
    - [alignment.py](tracetrack/entities/alignment.py) Contains Alignment class
    - [aligned.py](tracetrack/entities/aligned.py) Contains AlignedTrace and AlignedReference classes
    - [errors.py](tracetrack/entities/errors.py) Contains exception classes
    - [record.py](tracetrack/entities/record.py) Contains classes Feature, Record, TraceSeqRecord and find_quality_ends function
    - [reference_db.py](tracetrack/entities/reference_db.py) Contains ReferenceDb class and methods from reading reference sequences from files
    - [sequence_position.py](tracetrack/entities/sequence_position.py) Contains SequencePosition class
  - [templates](tracetrack/templates) 
    - [index.html](tracetrack/templates/index.html)  Source for main page of application, contains functionality for uploading files 
    - [layout.html](tracetrack/templates/layout.html) General layout of page
    - [help.html](tracetrack/templates/help.html) Help page
    - [loading.html](tracetrack/templates/loading.html) Page displayed while loading results
    - [results.html](tracetrack/templates/results.html) Source for page for displaying results table with alignments
    - [settings.html](tracetrack/templates/settings.html) Source of page for choosing alignment settings, such as quality trimming
    - [welcome.htlm](tracetrack/templates/welcome.html) Source for welcome page
  - [static](static)
    - [stylesheet.css](tracetrack/static/stylesheet.css) Defines style classes for use in the templates
    - [trace_viewer.js](tracetrack/static/trace_viewer.js) Script for showing trace file chromatograms
    - [table_display.js](tracetrack/static/table_display.js) Javascript functions for setting up results table
- [data](data)
   - [example](data/example) Contains example trace files and reference sequence tables
- [tests](tests) 
  - [test_alignment.py](tests/test_alignment.py) Tests for alignment module
  - [test_loading_data.py](tests/test_loading_data.py) Tests for loading_data module
- [resources](resources) Data for tests
- [environment.yml](environment.yml)
- [docker-compose.yml](docker-compose.yml) Configuration of task queue container
- [setup.py](setup.py) Script for tracetrack package installation
- [jenkins](jenkins)


## User Guide

For detailed instructions while use, see the "Help" tab in the application.

### Uploading files
Trace files can be uploaded in .ab1 or .zip format. If matching to reference by ID is desired, the trace file name has to contain the corresponding reference ID. Otherwise the file name can be arbitrary.  
The reference file can be an Excel spreadsheet (.xlsx file), a .csv file, a fasta file or a GenBank file. In case of csv and xlsx, the first row contains a header. Recommended column names are "ID" and "Sequence". A warning will be displayed otherwise, but the functionality is not affected as long as there is a header. The first column contains reference IDs, which can be used to match corresponding trace files. The second column contains the sequences. The reference IDs should be unique for each row. For reference sequences in fasta format, the ID in the header (after ">") is used as sequence ID. Similarly, the ID in a GenBank file is also used.

### Settings

#### Trace File Score Threshold
Each position of the sequence contained in a trace file is assigned a score of confidence. A slider enables the user to choose a score threshold. Any positions in the trace file with a lower score are discarded, the coverage in that part of the alignment will then be lower.

#### End Trimming Score Threshold
As the quality towards the ends of trace files tends to decrease, both ends of sequences are discarded until three bases in a row pass the trimming threshold. This can also be set by the user and should be higher than the score threshold, otherwise it has no effect.

#### Threshold for calling mixed peaks (*f*)
The application detects positions, where more than one trace peak is present, resulting in one primary and one or more secondary base calls.
This threshold specifies the minimal fraction of the area of the primary peak that another peak has to attain in order to be considered a secondary peak.

### Matching to Reference
Trace files are matched to reference sequences automatically. TraceTrack first checks if the reference ID is contained in the trace file name, and if no match is found this way, then the reference is chosen by best alignment score after aligning the trace file sequence to all the uploaded reference sequences. The read directionality is determined automatically. Both the reference and the directionality can be changed by the user manually.

### Data from Multiple Groups
If you wish to upload trace files that share the same reference ID, but come from a different group and you wish them to be aligned separately, upload one set of trace files using the available Choose Files button and another button will appear. Use this to upload files from other groups. You can repeat this step for multiple groups. Within each of the uploaded groups the trace files will be sorted by reference IDs and an alignment will be created for each ID in each group.

### Alignments
A multiple sequence alignment is created for each reference and its associated trace sequences. Trace files from each group are aligned separately. Substitution mutations, deletion mutations and insertion mutations are accepted when all trace files agree on them. In detail, the consensus sequence is produced by evaluating each alignment position separately using the following rules:

- If no reads cover the given position, the reference is kept, coverage is 0.</li>
- If at least one read matches the reference at the position, the reference nucleotide is kept and the coverage corresponds to the number of reads with this nucleotide at the given position.</li>
- If all the reads contain the same nucleotide, which differs from the reference, it is accepted as a mutation and the coverage is the number of reads at that position.
  - This also applies to the situation where there is only one read covering the given position.</li>
  - Gaps in trace files (deletion mutations) are handled same as substitutions - consensus between all trace files is required.</li>
- If all the aligned reads contain a gap in the reference (insertion mutation) and all of them contain the same base at given position, the insertion is accepted.</li>
- If some trace sequences contain a gap in the reference, an unknown insertion is displayed: "?".</li>
- If the reads differ from the reference but also from each other, the reference nucleotide is kept and coverage is 0.</li>

Every reference must contain at least one coding region (CDS feature). For reference sequences in GenBank format, features are extracted from the file. When reference sequences are uploaded in a different format, the whole sequence is considered coding. Coding and non-coding regions are displayed differently:
- Non-coding regions are displayed with lower opacity. CDS features have full opacity. Mutations are highlighted the same way.</li>
- Codon highlighting starts from the beginning of each CDS. Codons are not highlighted in the rest of the sequence.</li>
- If an insertion or deletion causes frameshift (shift by a number of positions not divisible by 3), the shifted region is highlighted in yellow and mutation types within the region are disregarded.</li>
- Frameshift is only displayed in CDS regions and is calculated from the beginning of the corresponding CDS.</li>

| Read      | base  |       |       |       |
|-----------|-------|-------|-------|-------|
| Reference | C     | C     | C     | C     |
| Read 1    | A     | A     | A     | A     |
| Read 2    | A     |       | A     | C     |
| Read 3    | A     |       | T     | A     |
| Result    | **A** | **A** | **C** | **C** |
| Coverage | 3     | 1     | 0     | 0     |

### Mixed peak detection
Each trace file is searched for heterogeneous positions. Signal to noise ratio is calculated for all positions as the ratio of the area of the main peak and the sum of areas of all other peaks. Areas with a low StN ratio are excluded from further steps, as these are noisy areas with low quality. Then an approximation of the area under curve fo each of the four traces is calculated at a given position, as well as heights of individual peaks. When a different base than the main called one has both peak area and height greater than a threshold (set to 15% of the main peak) and it's trace is concave around the center of the main peak, the position is marked as mixed. In the alignment, a position is marked as mixed if the secondary peak is detected in two or more trace files, or if it is detected in the only present trace file.

### Alignment Properties Table
The following properties of each alignment are displayed in the results table:
- Reference ID</li>
- % Coverage: What percentage of the reference sequence has coverage at least 1.</li>
- % Identity: What percentage of the reference nucleotides are exactly the same as the resulting ones.</li>
- Mutation Type: Which most severe type of mutation is present (none, silent, missense, nonsense).</li>
- Mutation counts: The numbers of each of the mutation types present.</li>
- Number of reads: How many trace files were associated with this reference.</li>
- File Names: which trace files were associated with this reference.</li>

The table can be ordered by different columns when the header is clicked.

### Displaying the Alignments
Each alignment can be displayed by clicking its Reference ID. The first character in the alignment (an asterisk) represents the Kozak sequence (GCCACC). If it is black, the intact sequence is present in at least one of the reads. If it is red, the Kozak sequence is either not covered by any read or mutated. The last character is an asterisk representing a stop codon. The same color coding applies as with the Kozak sequence. In the alignment itself, the background color represents coverage. Style of the letters represents mutations. Hovering your cursor over any nucleotide will highlight the codon with black and red underlining. If gaps are present in the alignment, black underlining represents the original reading frame and red is the frame after including gaps. A tooltip shows the reference nucleotide in the first row and individual reads in the following row.

### Download
The results can be downloaded as an Excel spreadsheet using the Download all button. The first sheet of the file will contain the same information as the results table. The reference IDs are clickable and take you to sheets containing the given alignments. The alignment sheet (named like the trace file, with "Seq" appended) contains a row for both the reference and consensus sequence, as well as each trace sequence. Amino acid translations are provided for each of these. The following sheet (marked by "Mut" appended to the sheet name)lists all the mismatched positions and regions of zero coverage. Position numbers are clickable and contain links to the corresponding position in the alignment.

If desired, each alignment can be downloaded as a separate spreadsheet using the button in the results table row.
