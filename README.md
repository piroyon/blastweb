# blastweb

Local BLAST search via Web UI and REST API using Flask, inspired by NCBI's legacy wwwblast.

- Search form with database / task / option choices filtered by program
- Multi-FASTA queries and FASTA file upload
- Searches run as background jobs; each result has its own URL
- wwwblast-style results: graphical overview, hit table and pairwise alignments
- Download results as pairwise text, TSV, CSV, XML, JSON or ASN.1 archive

---

### 🐍 Requirements

- Python 3.8 or later
- NCBI BLAST+（ include `makeblastdb`, `blastdbcmd`）
  - BLAST+ can be downloaded from here https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html


### 📦 Installation

Prepare venv: (Highly recommended)
```bash
python -m venv my_blastweb
source my_blastweb/bin/activate
```
from github: - If you choose this installation method, the ```blastweb``` command must be run in the ```<blastweb top dir>```.
```
git clone https://github.com/piroyon/blastweb.git
cd blastweb
pip install -e .
```
or, use pip:
```
pip install blastweb
```

### ⚙️ Configuration
Create a config file 'blast.yaml' in the current directory
```bash
blastweb init
```
Edit blast.yaml
```
blast_path: /usr/local/ncbi/blast+/bin
blast_db: /mnt/data/blastdb
default_extra_args: "-soft_masking true -word_size 11"
url_prefix: "/blastplus"
```

* ```blast_path```: Path to the directory where ```blastn```, ```blastp```, etc. are installed
* ```blast_db```: Path to the directory containing your BLAST databases
* ```default_extra_args```: Common BLAST options passed to all queries
* ```url_prefix```: If you plan to host this application behind a reverse proxy at a subpath
* ```timeout```: (optional) Max seconds a BLAST search may run. Default 600
* ```job_dir```: (optional) Where search jobs are stored. Default `<system temp dir>/blastweb_jobs`. Set this for production, since temp dirs may be cleaned
* ```job_retention_days```: (optional) Days to keep job results. Default 7
* ```max_running_jobs```: (optional) Max searches running at the same time. Default 4
* ```max_upload_mb```: (optional) Max size of an uploaded FASTA file. Default 10

Only common search options (e.g. `-word_size`, `-strand`, `-dust`, `-task`) are accepted in "Additional Options".
Options touching files, output format or remote search (e.g. `-out`, `-outfmt`, `-remote`) are rejected.
`default_extra_args` is set by the admin and is not restricted.

### 🚀 Start server for development or very local use

```bash
cd <blastweb_topdir>
blastweb runserver --port 5000 --config /path/to/blast.yaml
```
Then open http://localhost:5000 in your browser.

### 🧪 REST API
#### Synchronous search
POST to ```/api/blast``` with JSON:
```json
{
  "sequence": "ATGGCGTACGTAGC",
  "program": "blastn",
  "database": "mydb",
  "extra_args": "-word_size 11"
}
```
Command:
```bash
curl -X POST http://localhost:5000/api/blast \
  -H "Content-Type: application/json" \
  -d '{"sequence": "ATGGCGTACGTAGC", "program": "blastn", "database": "mydb"}'
```
Response (tsv lines split into array of columns):
```json
{
  "results": [["query1", "subject1", "98.7", "123", ...]]
}
```
Optional keys: `evalue`, `max_target_seqs`, `matrix`, `task`, `word_size`, `filter` (`yes`/`no`), `strand` (`plus`/`minus`).
If `extra_args` is omitted, `default_extra_args` is used.

#### Background jobs
For long searches, POST the same JSON to ```/api/jobs```. It returns at once with a job ID:
```json
{"job_id": "4df9...", "status_url": "http://localhost:5000/api/jobs/4df9...", "result_url": "http://localhost:5000/job/4df9..."}
```
Poll `status_url` until `state` is `done` (or `error`). When done, the response includes the BLAST JSON output (outfmt 15) as `blast_output`.
`result_url` is the result page in the browser.

### ✅ Running tests
Tests need BLAST+ on your PATH.
```bash
pip install -e ".[test]"
pytest
```

### 🗃️ Custom Databases
Prepare database for BLAST+ with `-hash_index` & `-parse_seqids` options.
```bash
makeblastdb -in mydb.fa -dbtype [nucl|prot] -hash_index -parse_seqids
```
Put your databases (e.g. ```mydb.fa.nin```, ```mydb.fa.nsq```, etc.) into the directory specified by ```blast_db```.
These will be auto-listed in the form as options, with their titles and types (via `blastdbcmd -list`). Alias databases (`.nal`/`.pal`) are listed too.

### 🔧 Production Deployment with Gunicorn and Nginx
For production environments, it is strongly recommended to run blastweb behind a WSGI server such as gunicorn and use a reverse proxy like nginx to handle TLS, compression, and static file delivery.
#### Example setup:
```css
[Client] → [nginx :443/80] → [gunicorn :5000] → blastweb.app
```
#### Step 1: Install gunicorn
```bash
pip install gunicorn
```
#### Step 2: Start blastweb with gunicorn
````bash
BLASTWEB_CONFIG=/path/to/blast.yaml gunicorn blastweb.wsgi:app --bind 127.0.0.1:5000
````
#### Step 3: Configure nginx (external to this project)
a minimal nginx config:
```nginx
server {
    listen 80;
    server_name your.domain.com;

    location / {
        proxy_pass http://127.0.0.1:5000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
    }
}
```
* To enable HTTPS, consider using [Let's Encrypt](https://letsencrypt.org/).
* Nginx is not included in this repository — please configure it separately on your server.

### 🏷️ Releasing (for maintainers)
Releases are published to PyPI by GitHub Actions (`.github/workflows/release.yml`) when a version tag is pushed.
1. Bump `version` in `pyproject.toml` and commit
2. Tag and push (the tag must match the version):
```bash
git tag v0.2.0
git push origin main v0.2.0
```
The workflow builds the package, uploads it to TestPyPI, then PyPI, and creates a GitHub Release with auto-generated notes.
