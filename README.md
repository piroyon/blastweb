# blastweb

Local BLAST search via Web UI and REST API using Flask.

---

## 🐍 Requirements

- Python 3.8 over
- NCBI BLAST+（ include `makeblastdb`）
  - BLAST+ can be downloaded from here https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html


## 📦 Installation

```bash
python -m venv my_blastweb
source my_blastweb/bin/activate
git clone https://github.com/<yourname>/blastweb.git
cd blastweb
pip install -e .
```

## ⚙️ Configuration
```bash
cp blast.yaml.example blast.yaml
```
Edit blast.yaml
```
blast_path: /usr/local/ncbi/blast+/bin
blast_db: /mnt/data/blastdb
default_extra_args: "-soft_masking true -word_size 11"
```

* ```blast_path```: Path to the directory where ```blastn```, ```blastp```, etc. are installed
* ```blast_db```: Path to the directory containing your BLAST databases
* ```default_extra_args```: Common BLAST options passed to all queries

## 🚀 Usage

Start the web server:
```bash
blastweb --port 5000
```
Then open http://localhost:5000 in your browser.

## 🧪 REST API
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

## 🗃️ Custom Databases
Put your databases (e.g. ```mydb.nin```, ```mydb.nsq```, etc.) into the directory specified by ```blast_db```.
These will be auto-listed in the form as options.
