# -*- coding: utf-8 -*-
from flask import Flask

app = Flask(__name__)

@app.route('/')
def index():
    return 'Hello Heroku_Flask'

if __name__ == '__main__':
    app.run()
