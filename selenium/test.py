from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select


def test():
	driver = webdriver.Chrome()
	driver.get("https://www.whois.com/whois")
 
	search_box = driver.find_element_by_css_selector('#query')  
	search_box.send_keys('www.google.com')               
	
	submit = driver.find_element_by_css_selector('#page-wrapper > div.whois-masthead > form > button')
	submit.click()
	driver.quit()


test()


