import argparse
import smtplib, ssl
from email.message import EmailMessage


def send_alert(args):
    """Send an email alert with the credentials passed in as args

    :param args: the command line args from argparse
    """

    message = EmailMessage()
    message.set_content(args.body)
    message['Subject'] = args.subject
    message['From'] = args.sender
    message['To'] = args.recipient

    context = ssl.create_default_context()
    with smtplib.SMTP(args.host, args.port) as mailserver:
        mailserver.starttls(context=context)
        mailserver.login(args.username, args.password)
        mailserver.send_message(message)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--host')
    parser.add_argument('--port')
    parser.add_argument('--username')
    parser.add_argument('--password')
    parser.add_argument('--sender')
    parser.add_argument('--recipient')
    parser.add_argument('--subject')
    parser.add_argument('--body')
    args = parser.parse_args()
    send_alert(args)