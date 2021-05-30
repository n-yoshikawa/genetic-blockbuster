# -*- coding: utf-8 -*-
import json

from flask import Flask, render_template, jsonify, request
from rdkit import Chem
from rdkit import rdBase

import crossover as co
import mutate as mu

# app.py
app = Flask(__name__)

lsmi = 'N#CC1=C(SCC(=O)Nc2cccc(Cl)c2)N=C([O-])[C@H](C#N)C12CCCCC2'
rsmi = 'O=[N+]([O-])c1c(Nc2cccc3ncccc23)ncnc1N1CCN(c2cccc(Cl)c2)CC1'
generation = 0

@app.route('/')
def index():
    return render_template("index.html")

@app.route('/start', methods=['POST'])
def start():
    return_json = {
        "leftSmiles": lsmi,
        "rightSmiles": rsmi,
        "generation": generation
    }
    return jsonify(values=json.dumps(return_json))

def mutate(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mutant = None
    while mutant is None:
        mutant = mu.mutate(mol, 1.0)
        try:
            Chem.SanitizeMol(mutant)
        except:
            continue
    return Chem.MolToSmiles(mutant)


@app.route('/mutate', methods=['POST'])
def show():
    global lsmi
    global rsmi
    global generation

    generation += 1
    smiles = request.json["smiles"]

    lsmi = mutate(smiles)
    rsmi = mutate(smiles)

    return_json = {
        "leftSmiles": lsmi,
        "rightSmiles": rsmi,
        "generation": generation
    }

    return jsonify(values=json.dumps(return_json))

if __name__ == '__main__':
    app.run(debug=True)
