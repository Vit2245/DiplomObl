import  requests
from bs4 import BeautifulSoup as BS
r=requests.get("https://sm.lsr.ru/spb/catalog/kirpich-gazobeton/kirpich/ryadovoy/?view=block")
html=BS(r.content,'html.parser')
for el in html.select(".b-product"):
    print(el)
    title=el.select('.b-title>a')
    print(title)