# -*- coding: utf-8 -*-
import copy
import nltk
import json

from flask import Flask, render_template, jsonify, request
import numpy as np

from rdkit import Chem
from rdkit import rdBase

import cfg_util
import zinc_grammar

rdBase.DisableLog('rdApp.error')
GCFG = zinc_grammar.GCFG
# app.py
app = Flask(__name__)

def CFGtoGene(prod_rules, max_len=-1):
    gene = []
    for r in prod_rules:
        lhs = GCFG.productions()[r].lhs()
        possible_rules = [idx for idx, rule in enumerate(GCFG.productions())
                          if rule.lhs() == lhs]
        gene.append(possible_rules.index(r))
    if max_len > 0:
        if len(gene) > max_len:
            gene = gene[:max_len]
        else:
            gene = gene + [np.random.randint(0, 256)
                           for _ in range(max_len-len(gene))]
    return gene


def GenetoCFG(gene):
    prod_rules = []
    stack = [GCFG.productions()[0].lhs()]
    for g in gene:
        try:
            lhs = stack.pop()
        except Exception:
            break
        possible_rules = [idx for idx, rule in enumerate(GCFG.productions())
                          if rule.lhs() == lhs]
        rule = possible_rules[g % len(possible_rules)]
        prod_rules.append(rule)
        rhs = filter(lambda a: (type(a) == nltk.grammar.Nonterminal)
                     and (str(a) != 'None'),
                     zinc_grammar.GCFG.productions()[rule].rhs())
        stack.extend(list(rhs)[::-1])
    return prod_rules


def mutation(gene):
    idx = np.random.choice(len(gene))
    gene_mutant = copy.deepcopy(gene)
    gene_mutant[idx] = np.random.randint(0, 256)
    return gene_mutant


def canonicalize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if smiles != '' and mol is not None and mol.GetNumAtoms() > 1:
        return Chem.MolToSmiles(mol)
    else:
        return smiles

lsmi = canonicalize('CC(=O)Nc1c2n(c3ccccc13)C[C@](C)(C(=O)NC1CCCCC1)N(C1CCCCC1)C2=O')
rsmi = canonicalize('O=[N+]([O-])c1c(Nc2cccc3ncccc23)ncnc1N1CCN(c2cccc(Cl)c2)CC1')
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

def get_mutant(smiles):
    new_smiles = ""
    gene = CFGtoGene(cfg_util.encode(smiles), max_len=300)
    while True:
        gene = CFGtoGene(cfg_util.encode(smiles), max_len=300)
        mutated_gene = mutation(gene)
        new_smiles = cfg_util.decode(GenetoCFG(mutated_gene))
        new_mol = Chem.MolFromSmiles(new_smiles)
        if len(new_smiles) < 15 or new_mol is None:
            continue
        try:
            Chem.SanitizeMol(new_mol)
        except:
            print("Sanitize failed:", new_smiles)
            continue
        new_smiles = Chem.MolToSmiles(new_mol)
        if new_smiles != smiles:
            break
    return new_smiles

@app.route('/mutate', methods=['POST'])
def show():
    global lsmi
    global rsmi
    global generation

    print(int(request.json["generation"]))
    if int(request.json["generation"]) > generation - 5:
        generation += 1
        smiles = canonicalize(request.json["smiles"])
        lsmi = get_mutant(smiles)
        rsmi = get_mutant(smiles)
        print('{},{},{}'.format(generation, lsmi, rsmi))

    return_json = {
        "leftSmiles": lsmi,
        "rightSmiles": rsmi,
        "generation": generation
    }

    return jsonify(values=json.dumps(return_json))


if __name__ == '__main__':
    app.run(debug=True)
