# readme
 In this directory code snippets are organized in separate markdown files based on software, task, etc.

 The [Wiki for this repo](https://github.com/RobertsLab/code/wiki) contains more detailed guides and instructions on things like:

 - [Jupyter Notebooks](https://github.com/RobertsLab/code/wiki/Jupyter-Notebook-Guide).
     - Jupyter Notebooks are an ideal means to document (and improve reproducibility) of computing commands.
 - [Super Computing with Mox (hyak)](https://github.com/RobertsLab/hyak_mox/wiki)
     - Mox is a high powered designed computer designed for seriously intensive computing
     - Due to the steep learning curve, it has its own dedicated wiki.
 - [Docker](https://github.com/RobertsLab/code/wiki/Docker-Guide)
     - Docker is a means to run a virtual Linux computer (great for Windows users!).
     - We have as Docker setup to install the most commonly used software used in the Roberts Lab, including Jupyter Notebook.

---

## File Descriptions

- **R.md** - Intended for R tips. Currently empty.

- **README.md** - This is the file you're currently reading.

- **bash.md** - Collection of commonly used bash commands for counting, finding & replacing, and moving/renaming files.

- **fasta.md** - FASTA file manipulation including counting sequences, converting FASTQ to FASTA, & converting FASTA to tab-delmited.

- **jupyter.md** - Ways to work with Jupyter Notebooks including embedding notebooks in websites, viewing notebooks on nbviwer.com, and using Jupyter Notebooks on remote computers.

- **misc.md** - Tips for using SSH, MD5 checksums, and gzip manipulation.

- `nf_core-base-srlab-500GB-node.config`: Nextflow config file for use with Roberts Lab 500GB Mox computing node.

- **[remote_connections.md](https://github.com/RobertsLab/code/blob/master/remote_connections.md)** - Ways to connect to computers remotely (via SSH), how to create/use SSH keys, create tunnels, and tmux to keep jobs running after closing SSH sessions.

- **snippets.cson**: Atom text editor snippets file which allows for autocompletion of markdown-formatted text frequently used by Roberts Lab members when writing in markdown (e.g. in notebooks). This file can be copied to `~/.atom/` on your computer for use and will overwrite the default snippets file.

- **sqlshare.md** - Tips for using SQLshare, but most commands are compatible with MySQL and/or SQLite. Lots of tips on joining tables and counting values.

- **vscode-snippets.code-snippets**: Visual Studio Code global snippets file which allows for autocompletion of markdown-formatted text frequently used by Roberts Lab members when writing in markdown (e.g. in notebooks). This file can be copied to `~/.config/Code/User/snippets/vscode-snippets.code-snippets` on your computer for use and will overwrite the default snippets file.
