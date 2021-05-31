# -*- coding: utf-8 -*-
import copy
import json
import random

from flask import Flask, render_template, jsonify, request, session
from rdkit import Chem
from rdkit import rdBase

import crossover as co
import mutate as mu

# app.py
app = Flask(__name__)
app.secret_key = 'hogefuga'

def choose_parents(population):
    lsmi = random.choice(population)
    population2 = copy.deepcopy(population)
    population2.remove(lsmi)
    rsmi = random.choice(population2)
    return lsmi, rsmi

@app.route('/')
def index():
    session['generation'] = 0
    session['population'] = ['CCO' for _ in range(5)] + ['c1ccccc1' for _ in range(5)]
    return render_template("index.html")

@app.route('/start', methods=['POST'])
def start():
    population = session['population']
    generation = session['generation']
    lsmi, rsmi = choose_parents(population)
    return_json = {
        "leftSmiles": lsmi,
        "rightSmiles": rsmi,
        "generation": generation
    }
    return jsonify(values=json.dumps(return_json))


@app.route('/reset', methods=['POST'])
def reset():
    session['generation'] = 0
    session['population'] = ['CCO' for _ in range(5)] + ['c1ccccc1' for _ in range(5)]
    population = session['population']
    generation = session['generation']
    lsmi, rsmi = choose_parents(population)
    return_json = {
        "leftSmiles": lsmi,
        "rightSmiles": rsmi,
        "generation": generation
    }
    return jsonify(values=json.dumps(return_json))


def reproduce(winner, loser):
    mol1 = Chem.MolFromSmiles(winner)
    mol2 = Chem.MolFromSmiles(loser)
    for _ in range(10):
        child = co.crossover(mol1, mol2)
        if child != None:
            child = mu.mutate(child, 0.5)
        try:
            Chem.SanitizeMol(child)
        except:
            continue
        smiles_child = Chem.MolToSmiles(child)
        if smiles_child == winner or smiles_child == loser:
            continue
        return smiles_child
    for _ in range(10):
        mutant = mu.mutate(mol1, 1.0)
        try:
            Chem.SanitizeMol(mutant)
        except:
            continue
        smiles_mutant = Chem.MolToSmiles(mutant)
        if smiles_mutant == winner or smiles_mutant == loser:
            continue
        return smiles_mutant
    return None


@app.route('/mutate', methods=['POST'])
def show():
    generation = session['generation']
    population = session['population']

    winner = request.json["winner"]
    loser = request.json["loser"]

    if winner in population and loser in population:
        child = reproduce(winner, loser)
        if child is not None and child not in population:
            generation += 1
            population.remove(loser)
            population.append(child)
            session['population'] = population
            session['generation'] = generation
    lsmi, rsmi = choose_parents(population)

    return_json = {
        "leftSmiles": lsmi,
        "rightSmiles": rsmi,
        "generation": generation
    }

    return jsonify(values=json.dumps(return_json))

if __name__ == '__main__':
    app.run(debug=True)
