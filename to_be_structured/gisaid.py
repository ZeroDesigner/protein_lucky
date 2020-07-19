import time
from selenium import webdriver

driver = webdriver.Chrome('path/to/chromedriver')  # Optional argument, if not specified will search path.
driver.get('https://platform.gisaid.org/epi3/frontend#2b8eee');
driver.refresh() #刷新页面
#driver.maximize_window() 
#填充用户名 密码 验证码
driver.find_element_by_id("elogin").send_keys('usename')
driver.find_element_by_id("epassword").send_keys('passwd')
driver.find_element_by_class_name("form_button_submit").click()
time.sleep(5)
driver.find_element_by_partial_link_text("EpiCoV™").click()
time.sleep(5)
driver.find_elements_by_class_name("sys-actionbar-action")[1].click()

page_num=67

for i in range(67):
    a=time.time()
    sim_num=len(driver.find_elements_by_class_name('yui-dt-rec'))
    for x in range(sim_num):
        print(x)
        driver.find_elements_by_class_name('yui-dt-rec')[x].click()
        time.sleep(3)
        driver.switch_to.frame(0)
        time.sleep(5)
        #meta下载
        driver.find_elements_by_class_name("sys-form-button-icon")[1].click()
        time.sleep(5)
        #fasta下载
        driver.find_elements_by_class_name("sys-form-button-icon")[2].click()
        time.sleep(10)
        driver.find_elements_by_class_name("sys-form-button-icon")[0].click()
        time.sleep(5)
    driver.find_element_by_class_name('yui-pg-next').click()
    time.sleep(5)
    b=time.time()
    print(b-a)
