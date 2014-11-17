#!/usr/bin/python
'''
Sending informations about result etc with Gmail
Using mail gmalert67@gmail.com
with password igbmcalert
'''

import smtplib
from email.MIMEMultipart import MIMEMultipart
from email.MIMEBase import MIMEBase
from email.MIMEText import MIMEText
from email import Encoders
import os

gmail_user = "gmalert67@gmail.com"
gmail_pwd = "igbmcalert"

def mail(to, subject, text= "", attach= None):
    
   # print "text ",text
   # print "attach ", attach
   msg = MIMEMultipart()

   msg['From'] = gmail_user
   msg['To'] = to
   msg['Subject'] = subject
   
   if text!='' :
       msg.attach(MIMEText(text))
   
   if attach is not None :
       
       part = MIMEBase('application', 'octet-stream')
       part.set_payload(open(attach, 'rb').read())
       Encoders.encode_base64(part)
       part.add_header('Content-Disposition',
               'attachment; filename="%s"' % os.path.basename(attach))
       msg.attach(part)

   mailServer = smtplib.SMTP("smtp.gmail.com", 587)
   mailServer.ehlo()
   mailServer.starttls()
   mailServer.ehlo()
   mailServer.login(gmail_user, gmail_pwd)
   mailServer.sendmail(gmail_user, to, msg.as_string())
   # Should be mailServer.quit(), but that crashes...
   mailServer.close()
if __name__ == '__main__': # Example
    mail("lionel.chiron@gmail.com","Hello from python!","This is a email sent with python","/Users/chiron/kurtosis_18.png")
