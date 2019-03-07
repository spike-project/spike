#!/usr/bin/env python 
# encoding: utf-8

#!/usr/bin/python
'''
Sending informations about result etc with Gmail
recipient (to) and attached document (attach) can be a list
'''

from __future__ import print_function
import smtplib
import os
import sys
if sys.version_info[0] < 3:
  from email.MIMEMultipart import MIMEMultipart
  from email.MIMEBase import MIMEBase
  from email.MIMEText import MIMEText
  from email import Encoders
else:
  from email import encoders as Encoders
  from email.mime.base import MIMEBase
  from email.mime.text import MIMEText
  from email.mime.multipart import MIMEMultipart

gmail_user = "gmalert67@gmail.com"
gmail_pwd = "IGBMCAlert67"

def add_to_msg(msg, f):
    part = MIMEBase('application', 'octet-stream')
    part.set_payload(open(f, 'rb').read())
    Encoders.encode_base64(part)
    part.add_header('Content-Disposition',
               'attachment; filename="%s"' % os.path.basename(f))
    msg.attach(part)
    return msg

def mail(to, subject, text= "", attach= None):
    
   # print "text ",text
   # print "attach ", attach
   msg = MIMEMultipart()
   msg['From'] = gmail_user
   if type(to) != list: # 
       msg['To'] = to
   else:                   # list of mails
       msg['To'] = ", ".join(to)
   msg['Subject'] = subject
  
   if text!='' :
      msg.attach(MIMEText(text))

   if attach is not None :
       #get all the attachments
       if type(attach) == list: # file exists, it is probably a list of files
           for f in attach:
               msg = add_to_msg(msg, f)
       else:
           msg = add_to_msg(msg, attach)

   mailServer = smtplib.SMTP("smtp.gmail.com", 587)
   mailServer.ehlo()
   mailServer.starttls()
   mailServer.ehlo()
   mailServer.login(gmail_user, gmail_pwd)
   mailServer.sendmail(gmail_user, to, msg.as_string())
   # Should be mailServer.quit(), but that crashes...
   mailServer.close()
if __name__ == '__main__': # Example
    mail("lionel.chiron@gmail.com","Hello from python!","This is a email sent with python","/Users/chiron/Pictures/kotok.jpg")
