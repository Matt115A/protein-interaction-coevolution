# utils/debug_entrez.py

import sys
import time
import logging
from Bio import Entrez

def debug_entrez(email=None, test_term="XP_007940319[Accession]", retries=3, delay=2):
    """
    Perform a quick Entrez esearch to verify the NCBI server is up.
    If it fails after `retries`, exits the program with an explanatory error.
    
    :param email: your Entrez.email; if provided, will be set on Entrez.
    :param test_term: an accession or term to query; default is the one you saw fail.
    :param retries: number of attempts before giving up.
    :param delay: base seconds to wait between attempts (will multiply by attempt number).
    """
    if email:
        Entrez.email = email

    for attempt in range(1, retries + 1):
        try:
            handle = Entrez.esearch(db="protein", term=test_term)
            _ = Entrez.read(handle)
            handle.close()
            logging.info(f"Entrez debug: successful esearch for '{test_term}'.")
            return True
        except Exception as e:
            logging.warning(f"Entrez debug attempt {attempt}/{retries} failed: {e}")
            if attempt < retries:
                time.sleep(delay * attempt)
            else:
                sys.stderr.write(
                    "\nERROR: Unable to reach NCBI Entrez server after "
                    f"{retries} attempts.\n"
                    "This indicates a problem on the server side. Please try again later.\n\n"
                )
                sys.exit(1)
