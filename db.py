
from sqlalchemy import create_engine
import pandas as pd
import glob
engine = create_engine('sqlite:///Inputs.db', echo=False)

csvfileList = glob.glob('xls/*.csv')

for filename in csvfileList:
    dfname = filename.split('\\')[1].split('.')[0]
    df = pd.read_csv(filename)
    df.to_sql(dfname, con=engine, if_exists='replace')


excelFileList = glob.glob('xls/*.xlsx')
available_sheets \
        = ['default_constants', 'country_constants', 'default_programs', 'country_programs', 'bcg_2014', 'bcg_2015',
           'bcg_2016', 'rate_birth_2014', 'rate_birth_2015', 'life_expectancy_2014', 'life_expectancy_2015',
           'notifications_2014', 'notifications_2015', 'notifications_2016', 'outcomes_2013', 'outcomes_2015',
           'mdr_2014', 'mdr_2015', 'mdr_2016', 'laboratories_2014', 'laboratories_2015', 'laboratories_2016',
           'strategy_2014', 'strategy_2015', 'strategy_2016', 'diabetes', 'gtb_2015', 'gtb_2016', 'latent_2016',
           'tb_hiv_2016', 'spending_inputs']

for filename in excelFileList:
    xls = pd.ExcelFile(filename)
    # Now you can list all sheets in the file

    if(len(xls.sheet_names) == 1):
      df_name = xls.sheet_names[0]
      print(df_name)
    else:
      numSheets = 0
      while numSheets < len(xls.sheet_names):
          #print(xls.sheet_names[numSheets])
          sheet_name = xls.sheet_names[numSheets]
          if sheet_name in available_sheets:
              #df_name = filename.split('\\')[1].split('.')[0] + '_' + sheet_name
              if sheet_name == 'rate_birth_2015' or sheet_name == 'life_expectancy_2015' :
                  df = pd.read_excel(filename, sheet_name=sheet_name, header = 3)
              else:
                  df = pd.read_excel(filename, sheet_name=sheet_name)
              print(sheet_name)
              df.to_sql(sheet_name, con=engine, if_exists='replace')
          numSheets = numSheets + 1

    def dbQuery(table_name, filter="", value="", column='*'):

        query ="Select "  + column + " from  " + table_name
        if filter != '' and value != '':
            query = "Select " + column + " from  " + table_name + "  Where " + filter + " =  \'" + value + '\''
        output_from_db = pd.read_sql_query(query, con=engine)
        return output_from_db


    res = dbQuery('bcg_2015', filter='Cname', value='Bhutan')
    print(res.values)
