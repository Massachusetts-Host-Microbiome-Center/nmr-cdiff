# NMR pre-processing
## Note
See [README.md](../README.md) in the base directory for dependencies before attempting to run this script.
## Usage
1. Navigate to the directory containing the raw NMR files from bruker using `cd`.
2. Run this script using `python /path/to/process.py`.
3. The script will prompt you for several fields. Simply follow the instructions.

Prior to running the script, it is recommended to place a text file `cfg.txt` in the raw NMR file directory, 
containing the reference ppm shifts and names of the expected chemicals, one per line, tab-separated. For example:

```
177.336 Proline
64.0059	Proline
48.8949	Proline
31.632	Proline
26.426	Proline
185.746	5-aminovalerate
42.1713	5-aminovalerate
39.4304	5-aminovalerate
29.2843	5-aminovalerate
25.2451	5-aminovalerate
```
[HMDB](https://hmdb.ca/) is a good source for reference chemical shifts.
