import requests
import socket
import os
import string
import random
import threading
import time
import psutil
import datetime
from bs4 import BeautifulSoup
from stem.control import Controller, Signal


# designed to create a text document containing specie and associated data
# this can then be read by a R script to create a table
class BacDive:

    def __init__(self):
        # base url address
        self.url = "https://bacdive.dsmz.de/strain/"
        # html content
        self.soup = None
        # all data found about microbe
        self.microbe_data = []
        # if next line contains value to desired data key
        self.next_contains = False
        # data key that found
        self.data_key = ""
        # data_value that found
        self.data_value = ""
        # which link currently on
        self.num = 0
        # get first page
        self.get_page()
        # file to save microbial data to
        self.visited_links_file_name = None

        # if using tor
        self.tor_address = '"C:/Users/micha/Desktop/Documents/Tor Browser/Browser/TorBrowser/Tor/tor.exe"'
        self.tor_session = None
        self.thread1 = None

    # get webpage
    def get_page(self):
        url = self.get_next_page(self.num)
        result = requests.session().get(url).text
        self.soup = BeautifulSoup(result, 'html.parser')
        self.num = self.num + 1

    # get desired specie and traits (have to manually get names first)
    def extract_info(self):
        # taxanomic data that want to retrieve
        taxonomic_trait = ["Domain", "Phylum", "Class", "Order",
                           "Family", "Genus"]
        # physiological data that want to retrieve
        physiology_trait = ["Gram stain", "Cell length", "Cell width",
                            "Cell shape", "Motility",
                            "Type of hemolysis", "Colony color",
                            "Colony size", "Colony shape",
                            "Incubation period",
                            "Ability of spore formation",
                            "Type of spore",
                            "Multicellular complex forming ability",
                            "Name of produced compound",
                            "Murein short key", "Murein types",
                            "Oxygen tolerance",
                            "Observation",
                            "Enzyme", "Nutrition type"]
        # data that want to retrieve that is found in a table format
        table_trait = ["met_antibiotica", "halophily",
                       "met_util", "met_production", "enzymes"]
        # growth conditions that want to retrieve
        growth_trait = ["Temperature range", "pH"]

        # first determine which specie
        specie = self.soup.select("body main div div div div div div div p")
        specie = BeautifulSoup(str(specie[3]), "lxml")
        specie = specie.span.get_text()
        self.microbe_data.append("|" + specie)
        print(specie)

        # get all rows that contain desired information
        rows = self.soup.find_all("td")
        for row in rows:
            row = BeautifulSoup(str(row), "lxml")
            row_text = row.td.get_text()

            # if find desired data key, next line will contain value
            if self.next_contains:
                # Append to list specie name, key, value
                self.data_value = row_text
                # add data
                self.microbe_data.append(self.format_data())
                self.next_contains = False
            else:
                # Search if contains key from lists
                if row_text in taxonomic_trait or row_text in physiology_trait or row_text in growth_trait:
                    self.data_key = row_text
                    self.next_contains = True

        # extract data in tables separately
        tables = self.soup.find_all("table")
        for table in tables:
            table = BeautifulSoup(str(table), "lxml")
            for trait in table_trait:
                wanted = table.find_all(attrs={"data-src-tbl": trait})
                for want in wanted:
                    self.data_key = trait
                    self.data_value = ""
                    child = want.parent.findChildren()
                    # Value(3) and whether trait expressed or not(5)
                    self.data_value = child[3].get_text() + child[5].get_text()
                    # add data
                    self.microbe_data.append(self.format_data())

    # Format data key and value when adding to file
    def format_data(self):
        data = ";" + self.data_key + ":" + self.data_value
        data = str.strip(data).replace('\\n', '').replace('\\t', '')
        return data

    # find hyperlink to next page
    # Each labeled by https://bacdive.dsmz.de/strain/<NUMBER>
    def get_next_page(self, num):
        return self.url + str(num)

    # recognize when there are no pages left
    def has_next_page(self):
        request = requests.get(self.url)
        if request.status_code < 400:
            return True
        else:
            return False

    # save data gathered to a file
    def setup_log(self, folder_name):
        # Track visited links
        visited_links_file_name = str(datetime.datetime.now().strftime("%Y-%m-%d %H-%M-%S")) + ".txt"
        self.visited_links_file_name = folder_name + "/" + visited_links_file_name
        open(self.visited_links_file_name, "x")

    def write_log(self):
        # Write hyperlinks to external final
        visit_file = open(self.visited_links_file_name, "a")
        for microbe in self.microbe_data:
            visit_file.write((str(microbe)))
        visit_file.close()
        self.microbe_data = []

    # get webpage through tor
    def tor_page(self):
        self.renew_ip()
        url = self.get_next_page(self.num)
        result = self.tor_session.get(url, timeout=2, headers=self.new_header()).text
        self.soup = BeautifulSoup(result, 'html.parser')
        self.num = self.num + 1

    # start tor socket
    def start_tor(self):

        # if tor not started, start
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        not_connected = True
        try:
            s.bind(('127.0.0.1', 9050))  # Try to open port
        except OSError:
            print("bound")
            not_connected = False
        s.close()
        if not_connected:
            print("connecting")
            # first connect to tor network
            self.thread1 = threading.Thread(target=self.launch_tor)
            self.thread1.start()

            # should check tat connected
            time.sleep(10)

        print("started tor")
        self.tor_session = requests.session()
        self.tor_session.proxies = dict()
        self.tor_session.proxies['http'] = 'socks5h://localhost:9050'
        self.tor_session.proxies['https'] = 'socks5h://localhost:9050'

        self.tor_session.cookies.clear()

    # open tor socket
    def launch_tor(self):
        os.system(self.tor_address)

    @staticmethod
    def renew_ip():
        """
        Obtain a new IP for tor port
        """
        with Controller.from_port(port=9051) as controller:
            controller.authenticate()
            controller.signal(Signal.NEWNYM)

    # get a new header for tor session
    @staticmethod
    def new_header():
        """
        New header for session
        :return: headers
        """
        headers = dict()
        letters = string.ascii_lowercase
        headers['User-agent'] = ''.join(random.choice(letters) for i in range(random.randint(1, 20)))
        return headers

    # end tor session
    def kill_tor(self):
        for proc in psutil.process_iter():
            if proc.name() == "tor.exe":
                proc.terminate()
        self.thread1.join()


# testing
visit_file = open("example.txt", "r")
html = visit_file.read()
soup = BeautifulSoup(html, "html.parser")
driver = BacDive()
driver.soup = soup
driver.setup_log("Microbial_Data")
driver.extract_info()
driver.write_log()

# driver = BacDive()
# driver.setup_log("Microbial_Data")
#
# driver.get_page()
# driver.extract_info()
# driver.write_log()

# total_count = 0  # for testing
# count = 0
# while driver.has_next_page() and total_count < 21:
#     driver.get_page()
#     driver.extract_info()
#     count = count + 1
#     total_count = total_count + 1
#     time.sleep(.1)
#     if count >= 20:
#         driver.write_log()


# if using tor
# driver.start_tor()
# driver.tor_page()
# driver.extract_info()
# driver.kill_tor()
