from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
import time


def test():
	driver = webdriver.Chrome()
	driver.get("https://www.whois.com/whois")
 
	search_box = driver.find_element_by_css_selector('#query')  
	search_box.send_keys('www.google.com')               
	
	submit = driver.find_element_by_css_selector('#page-wrapper > div.whois-masthead > form > button')
	submit.click()
	
	register = driver.find_element_by_css_selector('#page-wrapper > div > div.whois_main_column > div:nth-child(6) > div:nth-child(5) > div.df-value')
	print(f"Registered On: {register.text}")
	
	expire = driver.find_element_by_css_selector('#page-wrapper > div > div.whois_main_column > div:nth-child(6) > div:nth-child(5) > div.df-value')
	print(f"Expires On: {expire.text}")
	
	time.sleep(10)
	driver.quit()
test()


