'''
Created by Lionel Chiron  18/10/2013 
Copyright (c) 2013 __NMRTEC__. All rights reserved.
'''
import sys, os
import threading
import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email import Encoders
from time import sleep, localtime, strftime
from uuid import getnode as get_mac
import os.path as op

gmail_user = "gmalert67@gmail.com"
gmail_pwd = "igbmcalert"

def mail(to, subject, text = "", attach = None):
    
   print "text ",text
   print "attach ", attach
   msg = MIMEMultipart()
   #################
   msg['From'] = gmail_user
   msg['To'] = to
   msg['Subject'] = subject
   #############
   if text!='' :
       msg.attach(MIMEText(text))
   ###############
   if attach is not None :
       part = MIMEBase('application', 'octet-stream')
       part.set_payload(open(attach, 'rb').read())
       Encoders.encode_base64(part)
       part.add_header('Content-Disposition',
               'attachment; filename="%s"' % os.path.basename(attach))
       msg.attach(part)
   ##############
   mailServer = smtplib.SMTP("smtp.gmail.com", 587)
   mailServer.ehlo()
   mailServer.starttls()
   mailServer.ehlo()
   mailServer.login(gmail_user, gmail_pwd)
   mailServer.sendmail(gmail_user, to, msg.as_string())
   # Should be mailServer.quit(), but that crashes...
   mailServer.close()

class Logger(object):
    def __init__(self):
        try:
            applic = sys.modules['__main__'].__file__
        except Exception:
            print "no __file__"
            applic = ''
        
        self.terminal = sys.stderr
        date = self.datetime()
        mac = self.mac()+'_'
        app = op.splitext(op.basename(applic))[0]+'_'
        self.log_name = "log_" + app + mac + date + ".dat"
        self.log = open(self.log_name, "w")
        self.trig = False
    
    def send_mail(self):
        a = threading.Thread(None, self.mail_if_error, None, )
        a.start()
    
    def mac(self):
        mac_addr = str(get_mac())
        mac_dic = {'149885691548389': 'kartikeya'} #149885691548389
        if mac_addr in mac_dic :
            return mac_dic[mac_addr]
        else:
            return mac_addr
    
    def datetime(self):    
        return strftime('%Y-%m-%d-%H-%M', localtime()) 
    
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        if not self.trig :
            self.send_mail()
            self.trig = True

    def mail_if_error(self):
        sleep(1)
        self.log.close()
        mail("lionel.chiron@gmail.com", self.log_name, attach = self.log_name)
        os.remove(self.log_name)

if __name__ == '__main__':
    sys.stderr = Logger()
    f = open("fff.jpg")

    