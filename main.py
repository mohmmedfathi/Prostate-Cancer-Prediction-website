import streamlit as st
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
import os
from sklearn.metrics import accuracy_score
import pickle
from st_aggrid import AgGrid, GridOptionsBuilder
from st_aggrid.shared import GridUpdateMode
from fpdf import FPDF
import base64
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq
from streamlit_option_menu import option_menu
import altair as alt

def fun ():
    pickle_in = open("classifier.pkl", "rb")
    classifier = pickle.load(pickle_in)
    return classifier
def prediction(x):
     prd = fun()
     prediction = prd.predict(x)
     return prediction[0]

with st.sidebar:
    selected = option_menu(menu_title='Main Menu',
                           options=['Home','Dataset Analysis','Model Prediction','DNA Sequence','DNA Manipulation'],
                           icons=['house','book','boxes','flower2','stars'],
                           menu_icon='cast',
                           default_index=0
                          # orientation='horizontal'
                           )

def HomePage():
    st.snow()


    st.title("Bioinformatics Project")
    st.header("Prof.Sara Salah")
    #st.header("Dr.")
    st.header("TEAM NAME : ")
    st.subheader("  Mohammed Fathi Sayed Ahmed")
    st.subheader("  Mohammed Medhat Mohammed ")
    st.subheader("  Mahmoud AbdElGhany")
    st.subheader("  Mahmoud Medhat")
    st.subheader("  Mohammed Misaraa")
    st.subheader("  Omar Alaa")


def DNA_nucleotide_count(seq):
        d = dict([
            ('A', seq.count('A')),
            ('T', seq.count('T')),
            ('G', seq.count('G')),
            ('C', seq.count('C'))
        ])
        return d
def DNAtool():



    st.write("""
    # DNA Nucleotide Tool
    This app counts the nucleotide of query DNA!
    ***
    """)

    ######################
    # Input Text Box
    ######################

    # st.sidebar.header('Enter DNA sequence')
    st.header('Enter DNA sequence')

    sequence_input = "GAACACGTGGAGGCAAACAGGAAGGTGAAGAAGAACTTATCCTATCAGGACGGAAGGTCCTGTGCTCGGG\nATCTTCCAGACGTCGCGACTCTAAATTGCCCCCTCTGAGGTCAAGGAACACAAGATGGTTTTGGAAATGC\nTGAACCCGATACATTATAACATCACCAGCATCGTGCCTGAAGCCATGCCTGCTGCCACCATGCCAGTCCT"

    # sequence = st.sidebar.text_area("Sequence input", sequence_input, height=250)
    sequence = st.text_area("Sequence input", sequence_input, height=250)
    #sequence = sequence.splitlines()
    #sequence = sequence[1:]  # Skips the sequence name (first line)
    #sequence = ''.join(sequence)  # Concatenates list to string

    #st.write("""
    #***
    #""")

    ## Prints the input DNA sequence

    ## DNA nucleotide count
    st.header('OUTPUT (DNA Nucleotide Count)')

    ### 1. Print dictionary
    #st.subheader('1. Print dictionary')

    X = DNA_nucleotide_count(sequence)
    st.info(X)

    ### 2. Print text
    ss ='There are  ' + str(X['A']) + ' adenine (A) \n\n There are  ' + str(X['T']) + ' thymine (T) \n \nThere are  '+str(X['G']) + ' guanine (G) \n\n There are  ' + str(X['C']) + ' cytosine (C)'

    st.subheader('1. Print text')
    st.info(ss)


    ### 3. Display DataFrame
    st.subheader('2. Display DataFrame')
    df = pd.DataFrame.from_dict(X, orient='index')
    df = df.rename({0: 'count'}, axis='columns')
    df.reset_index(inplace=True)
    df = df.rename(columns={'index': 'nucleotide'})
    st.write(df)

    ### 4. Display Bar Chart using Altair
    st.subheader('3. Display Bar chart')
    p = alt.Chart(df).mark_bar().encode(
        x='nucleotide',
        y='count'
    )
    p = p.properties(
        width=alt.Step(80)  # controls width of bar.
    )
    st.write(p)


def remBase(str, ch):
    ans = ''
    for i in str:
        if i != ch:
            ans += i
    return ans
def DNA_manipulation():
       x_upload = st.file_uploader("upload file")
       str = ''
       if x_upload:
           t = ''
           for line in x_upload:
               t += line.decode('utf-8')
           str += t
           st.subheader("File Content")
           st.write(str)
           new_title = '<p style="font-family:sans-serif; color:Blue; font-size: 22px;">Complement DNA</p>'
           st.markdown(new_title, unsafe_allow_html=True)
           codDna = Seq(str)
           st.write(codDna.complement())

           new_title = '<p style="font-family:sans-serif; color:Blue; font-size: 22px;">Transcribe DNA</p>'
           st.markdown(new_title, unsafe_allow_html=True)

           st.write(codDna.reverse_complement().transcribe())

           new_title = '<p style="font-family:sans-serif; color:Blue; font-size: 22px;">Translation DNA</p>'
           st.markdown(new_title, unsafe_allow_html=True)

           messengerRna = Seq(str, IUPAC.unambiguous_rna)
           st.write("messengerRna = " + messengerRna)
           st.write("Translated messengerRna = " + messengerRna.translate())

           edit = st.checkbox("do you want to edit ?")
           if edit:
               changeBase = st.checkbox("1 - Change Specific Base")
               if changeBase:
                   st.warning("enter index of base :")
                   index = st.number_input("index")
                   if index >= 0 and index < len(str):
                       st.warning("enter base (A,T,C,G):")
                       ch = st.text_input("base ")
                       if ch:
                           if not (ch == 'A' or ch == 'T' or ch == 'C' or ch == 'G'):
                               st.error('enter base only (A,T,C,G)')
                               # str[index] = ch
                           else:
                               s = list(str)
                               s[int(index)] = ch
                               "".join(s)
                               st.write(s)
                   else:
                       st.error('enter valid index.')
               removeBase = st.checkbox("2 - remove Specific Base")
               if removeBase:
                   rd = st.radio("Select base to entirely from Sequence ",
                                 ('A', 'T', 'C', 'G'))
                   mutable_obj = MutableSeq(str, IUPAC.unambiguous_rna)
                   flag = 0
                   st.warning("before Removal : " + mutable_obj)
                   if rd == "A":
                       xd = remBase(mutable_obj, 'A')
                       flag = 1
                   elif rd == 'C':
                       xd = remBase(mutable_obj, 'C')
                       flag = 1
                   elif rd == 'T':
                       xd = remBase(mutable_obj, 'T')
                       flag = 1
                   elif rd == 'G':
                       xd = remBase(mutable_obj, 'G')
                       flag = 1

                   if flag == 1:
                       st.success("after Removal: " + xd)




def ModelPredictionPage():
    st.title("Prostate Cancer Prediction Model")

    paitentId = st.text_input("Patient ID")
    radius = st.text_input(label='enter radius')
    area = st.text_input(label='enter area')
    smoothness = st.text_input(label='enter smoothness')
    compactness = st.text_input(label='enter compactness')
    symmetry = st.text_input(label='enter symmetry')

    if st.button('Predict'):
        if (prediction([[radius, area, smoothness, compactness, symmetry]]) == 0):
            #nnp = '<p style = "font-family:sans-serif;color:Red;font-size:42px;">the result is : Malignant Tumor</p>'
            nnp = 'Malignant Tumor'
            st.error(nnp)

        else:
            #nnp = '<p style = "font-family:sans-serif;color:White;font-size:42px;">the result is : Benign tumor</p>'
            nnp = 'Benign tumor'
            st.warning(nnp)

    pdf = st.checkbox("Export as pdf")
    if pdf:

     def create_download_link(val, filename):
        b64 = base64.b64encode(val)  # val looks like b'...'
        return f'<a href="data:application/octet-stream;base64,{b64.decode()}" download="{filename}.pdf">Download file</a>'


     pdf = FPDF()
     pdf.add_page()
     pdf.set_font('Arial', 'B', 16)
     pdf.image('icon.png', x = 140, y = 30, w = 50,h = 37)
     pdf.cell(0, 10, 'Prdiction for Patient ID : '+ str(paitentId), 0, 1,align='C') #

     pdf.cell(0, 10, 'radius : ' + str(radius),0,1)#
     pdf.cell(0, 10, 'area : ' + str(area), 0, 1)#
     pdf.cell(0, 10, 'smoothness : '+ str(smoothness) , 0, 1)#
     pdf.cell(0, 10, 'compactness : '+ str(compactness), 0, 1)#
     pdf.cell(0, 10, 'symmetry : '+ str(symmetry) , 0, 1) #
     if (prediction([[radius, area, smoothness, compactness, symmetry]]) == 0):
             nnp = 'the result is : Malignant Tumor'
     else:
             nnp = 'the result is : Benign tumor'

     if nnp=='the result is : Malignant Tumor':
            pdf.set_text_color(255, 0, 0)
     else:
            pdf.set_text_color(0, 255, 255)
     pdf.set_font('Arial', 'B', 25)
     pdf.cell(0, 10, nnp, 0, 1,align='C')
        #

     html = create_download_link(pdf.output(dest="S").encode("latin-1"), "test_result")

     st.markdown(html, unsafe_allow_html=True)

     #st.balloons()

def excuteInteractive():
    df = pd.read_csv("Prostate_Cancer.csv")
    gd = GridOptionsBuilder.from_dataframe(df)
    gd.configure_pagination(enabled=True)
    gd.configure_default_column(editable=True, groupable=True)

    sel_mode = st.radio('Selection Type', options=['single', 'multiple'])

    gd.configure_selection(selection_mode=sel_mode, use_checkbox=True)
    gridoptions = gd.build()
    grid_table = AgGrid(df, gridOptions=gridoptions, update_mode=GridUpdateMode.SELECTION_CHANGED,
                        height=500,
                        allow_unsafe_jscode=True,
                        reload_data=True,
                        theme='fresh')
    sel_row = grid_table['selected_rows']
    st.subheader("output")
    st.write(sel_row)
def DatasetAnalysisPage():


      interactiveTable = st.checkbox("interactive table of Prostate_Cancer Dataset")
      if interactiveTable:
        excuteInteractive()
      else:
        data = pd.read_csv("Prostate_Cancer.csv")
        st.header("Data processing")

       # st.subheader("Data info")
        #st.write(data.info())

        st.subheader("Data Stastistics")
        st.write(data.describe())

        st.subheader("first 5 patient in Dataset")
        st.write(data.head())

        st.subheader("last 5 patient in Dataset")
        st.write(data.tail())
        st.header("Data Visualization")
        st.subheader("Correlation")
        data['diagnosis_result'].replace({'M': 0, 'B': 1}, inplace=True)
        corr_metrics = data.corr()
        st.write(corr_metrics.style.background_gradient())

        # st.subheader("Heat map")
        # corr = data.corr()
        # plt.figure(figsize=(15, 10))
        # sns.heatmap(corr, annot=True, cmap=plt.cm.Reds)
        # plt.show()

    #username = st.text_input("Username")
    #password = st.text_input("Password",type = 'password')
    #st.button("Login")
if selected == 'Home':
    HomePage()
elif selected == 'Model Prediction':
    ModelPredictionPage()
elif selected == "Dataset Analysis":
    DatasetAnalysisPage()
elif selected == "DNA Sequence":
    DNAtool()
elif selected == 'DNA Manipulation':
    DNA_manipulation()

# # import pandas as pd
# # import streamlit as st
# # import altair as alt
# # from PIL import Image
import streamlit as st

