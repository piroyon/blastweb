# blastweb

Local BLAST search via Web UI and REST API using Flask.

---

### 🐍 Requirements

- Python 3.8 over
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
### 🚀 Start server for development or very local use

```bash
cd <blastweb_topdir>
blastweb runserver --port 5000 --config /path/to/blast.yaml
```
Then open http://localhost:5000 in your browser.

### 🧪 REST API
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

### 🗃️ Custom Databases
Prepare database for BLAST+ with `-hash_index` & `-parse_seqids` options.
```bash
makeblastdb -in mydb.fa -dbtype [nucl|prot] -hash_index -parse_seqids
```
Put your databases (e.g. ```mydb.fa.nin```, ```mydb.fa.nsq```, etc.) into the directory specified by ```blast_db```.
These will be auto-listed in the form as options.

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
