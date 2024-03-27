VENV_CREATED=10
echo $VENV_CREATED
(if [[ $VENV_CREATED == 10 ]]; then
    echo "a"
    AVENVA=20
    echo $AVENVA
fi) && echo $AVENVA
