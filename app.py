from flask import Flask, render_template, request
from CGViewBuilder import CGViewBuilder
import os

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route("/upload", methods=['POST'])
def upload():
    dataset_file = request.files.get('dataset')
    config_file = request.files.get('config')
    
    config_file.save(config_file.filename) if config_file else None

    if dataset_file:
        dataset_file.save(dataset_file.filename)

        cgview_options = {
            'config': config_file.filename,
            'contigs': None
        }
        data = CGViewBuilder(dataset_file.filename, options=cgview_options)
        
        os.remove(dataset_file.filename)
        os.remove(config_file.filename) if config_file else None
        
        return render_template("result.html", data=data.to_json())
    else:
        return "No file uploaded."
            

if __name__ == "__main__":
    app.run(debug=True)