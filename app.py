"""Front-end for patent search application."""

import datetime
import string

from flask import Flask, render_template, request
from nltk.stem.porter import PorterStemmer
import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import linear_kernel
import spacy
from spacy.lang.en import English
from spacy.lang.en.stop_words import STOP_WORDS


STEMMER = PorterStemmer()

# Load the spaCy english library, PUNCTUATION, and stopwords
NLP = spacy.load("en_core_web_sm")
# add a shortlist of patent stop-words
NLP.Defaults.stop_words |= {
    "electrochemical",
    "electrochemistry",
    "use",
    "comprise",
    "invention",
    "composition",
    "form",
    "include",
    "active",
    "method",
    "compound",
    "group",
    "high",
    "present",
    "provide",
    "produce",
    "contain",
}

PUNCTUATION = string.punctuation
PARSER = English()


# Read the patent dataset
DF_USER_SEARCH = pd.read_csv("data/user_search.csv")
DF = pd.read_csv("data/export_patent_view_main_tokenized_v2.csv")

# Define the TFIDF vectorizer
VECTORIZER = TfidfVectorizer(ngram_range=(1, 3))

# Imported dataset contains the patent abstracts pre-processed (using spacy_tokenizer function below) into tokens as a column.  Use that column as the features.
DOCS = DF["tokens"]
BAG = VECTORIZER.fit_transform(DOCS)

APP = Flask(__name__)


def spacy_tokenizer(user_doc):
    """Process the user's input into tokens.
    
    Process user input into tokens that can be compared against
    the tokenized patent abstract data
    
    Args:
        user_doc (str): The string of user input to tokenize
    
    Returns:
        Tuple[pandas.Series, pandas.Series]: Tuple with two Series:
            1. Series of outgoing tokens
            2. Series of tokens to print
    """

    mytokens = []
    mytokens_out = []
    print_tokens = []
    mytokens = PARSER(user_doc)
    # Process text NO1: lowercase, remove spaces, and lemmatize
    mytokens = [
        word.lemma_.lower().strip() if word.lemma_ != "-PRON-" else word.lower_
        for word in mytokens
    ]
    # Process text NO2: remove stop words and PUNCTUATION
    mytokens = [
        word for word in mytokens if word not in STOP_WORDS and word not in PUNCTUATION
    ]

    mytokens = NLP(str(mytokens))
    for word in mytokens:
        # remove verbs, punctuation, conjunctions, adverbs, numbers, and undefined parts-of-speech
        mytoken = (
            (word.pos_ != "VERB")
            & (word.pos_ != "PUNCT")
            & (word.pos_ != "CCONJ")
            & (word.pos_ != "ADV")
            & (word.pos_ != "NUM")
            & (word.pos_ != "X")
        )
        if mytoken:
            print_tokens.append(word)
            # stem
            # word=STEMMER.stem(word.text)
            mytokens_out.append(word)

    mytokens_out = pd.Series(str(mytokens_out))
    return mytokens_out, print_tokens


# initial page front
@APP.route("/")
def home():
    return render_template("index.html.j2")


# page that returns after a user "submits"
@APP.route("/analyze", methods=["POST"])
def analyze():
    if request.method == "POST":
        # Bring in the text query and number of results wanted by the user
        rawtext = request.form["rawtext"]
        n = request.form["result_count"]
        n = int(n)
        received_text = rawtext

        # store user search for analysis
        DF_USER_SEARCH = pd.read_csv("data/user_search.csv")
        df_store = pd.DataFrame(
            {"timestamp": [datetime.datetime.now()], "user_input": [rawtext]}
        )
        DF_USER_SEARCH = pd.concat(
            [DF_USER_SEARCH, df_store], ignore_index=True, sort=True
        )
        DF_USER_SEARCH = DF_USER_SEARCH[["timestamp", "user_input"]]
        DF_USER_SEARCH.to_csv("data/user_search.csv")

        # convert input into tokens for the user to see
        text_tokens, print_tokens = spacy_tokenizer(received_text)
        # vectorizer user input
        user_input = pd.Series(text_tokens)
        # transform tokenized query into TFIDF
        test = VECTORIZER.transform(user_input)
        # test cosine similarities against all the documents
        cosine_similarities = linear_kernel(test, BAG).flatten()
        # return indices of most similar documents
        related_docs_indices = cosine_similarities.argsort()[: -(n + 1) : -1]
        # store similarity in the dataframe
        DF["similarity"] = cosine_similarities

        # extract the top n to be returned to the index.html
        out_df = DF.iloc[list(related_docs_indices)][
            [
                "title",
                "abstract",
                "patent_id",
                "similarity",
                "year",
                "subclass_id",
                "url",
            ]
        ]
        # reset the indeces for the returned data
        out_df = out_df.reset_index()
        # convert dataframe into dictionary
        out_df = out_df.to_dict("records")

    return render_template(
        "analyze.html.j2",
        received_text=rawtext,
        text_tokens=print_tokens[:],
        table=out_df,
    )


# https://www.websiteout.net/counter.php


if __name__ == "__main__":
    APP.run(debug=True, port=5957)
