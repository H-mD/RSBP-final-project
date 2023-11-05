from functools import wraps
from flask import Flask, render_template, request, redirect, session, url_for
from CGViewBuilder import CGViewBuilder
from flask_mysqldb import MySQL
import MySQLdb.cursors
import re
import os

app = Flask(__name__)

app.secret_key = 'hafidganteng'

app.config['MYSQL_HOST'] = 'localhost'
app.config['MYSQL_USER'] = 'hafid'
app.config['MYSQL_PASSWORD'] = '1234'
app.config['MYSQL_DB'] = 'rsbp2'

mysql = MySQL(app)

# Define the login_required decorator
def login_required(f):
    @wraps(f)
    def decorated_function(*args, **kwargs):
        if 'loggedin' in session:
            return f(*args, **kwargs)
        else:
            return redirect(url_for('login'))
    return decorated_function

@app.route('/')
@app.route('/login', methods=['GET', 'POST'])
def login():
    message = ''
    if request.method == 'POST' and 'email' in request.form and 'password' in request.form:
        email = request.form['email']
        password = request.form['password']
        cursor = mysql.connection.cursor(MySQLdb.cursors.DictCursor)
        cursor.execute('SELECT * FROM user WHERE (email = %s) AND password = %s', (email, password,))
        user = cursor.fetchone()
        if user:
            session['loggedin'] = True
            session['userid'] = user['userid']
            session['name'] = user['name']
            session['email'] = user['email']
            message = 'Logged in successfully!'
            return redirect(url_for('index'))
        else:
            message = 'Please enter the correct email / password!'
    return render_template('login.html', message=message)

@app.route('/logout')
@login_required
def logout():
    session.pop('loggedin', None)
    session.pop('userid', None)
    session.pop('email', None)
    return redirect(url_for('login'))

@app.route('/register', methods=['GET', 'POST'])
def register():
    message = ''
    if request.method == 'POST' and 'name' in request.form and 'password' in request.form and 'email' in request.form:
        userName = request.form['name']
        password = request.form['password']
        email = request.form['email']
        cursor = mysql.connection.cursor(MySQLdb.cursors.DictCursor)
        cursor.execute('SELECT * FROM user WHERE email = %s ', (email,))
        account = cursor.fetchone()
        if account:
            message = 'Account already exists!'
        elif not re.match(r'[^@]+@[^@]+\.[^@]+', email):
            message = 'Invalid email address!'
        elif not userName or not password or not email:
            message = 'Please fill out the form!'
        else:
            cursor.execute('INSERT INTO user VALUES (NULL, %s, %s, %s)', (userName, email, password,))
            mysql.connection.commit()
            message = 'You have successfully registered!'
            return render_template('login.html', message=message)
    elif request.method == 'POST':
        message = 'Please fill out the form!'
    return render_template('register.html', message=message)

@app.route('/index')
@login_required
def index():
    return render_template('index.html')

@app.route("/upload", methods=['POST'])
@login_required
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
