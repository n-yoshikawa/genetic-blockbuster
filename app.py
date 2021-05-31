# -*- coding: utf-8 -*-
import copy
import json
import random

from flask import Flask, render_template, jsonify, request
from rdkit import Chem
from rdkit import rdBase

import crossover as co
import mutate as mu

# app.py
app = Flask(__name__)

generation = 0
population = ['CCO' for _ in range(5)] + ['c1ccccc1' for _ in range(5)]

def choose_parents():
    lsmi = random.choice(population)
    population2 = copy.deepcopy(population)
    population2.remove(lsmi)
    rsmi = random.choice(population2)
    return lsmi, rsmi

@app.route('/')
def index():
    return render_template("index.html")

@app.route('/start', methods=['POST'])
def start():
    lsmi, rsmi = choose_parents()
    return_json = {
        "leftSmiles": lsmi,
        "rightSmiles": rsmi,
        "generation": generation
    }
    return jsonify(values=json.dumps(return_json))


def reproduce(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles1)
    for _ in range(100):
        child = co.crossover(mol1, mol2)
        if child != None:
            child = mu.mutate(child, 0.5)
        try:
            Chem.SanitizeMol(child)
        except:
            continue
        smiles_child = Chem.MolToSmiles(child)
        if smiles_child == smiles1 or smiles_child == smiles2:
            continue
        return smiles_child
    return None


@app.route('/mutate', methods=['POST'])
def show():
    global generation
    global population

    winner = request.json["winner"]
    loser = request.json["loser"]

    if winner in population and loser in population:
        child = reproduce(winner, loser)
        if child is not None and child not in population:
            generation += 1
            population.remove(loser)
            population.append(child)
    lsmi, rsmi = choose_parents()

    return_json = {
        "leftSmiles": lsmi,
        "rightSmiles": rsmi,
        "generation": generation
    }

    return jsonify(values=json.dumps(return_json))

if __name__ == '__main__':
    app.run(debug=True)
