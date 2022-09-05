# Scripts for generation of random barcodes

Each python script takes as only argument a config file containing all the parameters to:

- Generate barcodes (1_BarcodeGen_split.py)
- Simulate garbled reads from the list of barcodes (2_SimGen_split.py)
- Demultiplexing the reads (3_Barcodedecoding)
- Calculate the precision and recall after demultiplexing of the reads (4_BarcodeRecall.py)

To run these scripts a python module including Pythorch must be load

```
module load PyTorch/1.10.0-foss-2021a-CUDA-11.3.1
```
Also the library networkx must be installed

```
pip install networkx
```

Run the desired script like this:

```
python 1_BarcodeGen_split.py --config "config_file.yaml"
```

